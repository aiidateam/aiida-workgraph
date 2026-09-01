from __future__ import annotations

from enum import Enum
from typing import Any, Dict, Optional

from aiida_pythonjob.data.serializer import all_serializers
from aiida_pythonjob.utils import serialize_ports
from node_graph.serializer import SerializationAdapter
from node_graph.utils import resolve_tagged_values


def _flatten_enums(value: Any) -> Any:
    """Replace ``Enum`` members with their bare ``.value``, recursively.

    ``general_serializer`` (aiida-pythonjob) has no serializer for ``Enum``
    members: a raw ``Enum`` -- or one nested in a dict/list/tuple -- reaches
    it and the whole submission fails at submit time with an opaque
    serialization error. Collapsing members to their ``.value`` here keeps
    the payload JSON-serializable so ``serialize_ports`` succeeds. The
    ``isinstance`` check also matches wrapt proxies (node-graph's
    ``TaggedValue``) wrapping an enum member.

    Flattening is one-way. Whether a body then receives the ``Enum`` it
    declared is the installed ``node_graph``'s call: the read side
    (``coerce_inputs_from_spec``) rebuilds a socket's value only from the
    ``structured_type`` descriptor its spec records, and not every
    ``node_graph`` records one for ``Enum``. Without it a function task's
    body receives the bare value and a ``@task.graph`` body the stored
    ``orm.Str``; with it both receive the member, a ``@task.graph`` body
    behind a ``TaggedValue``. A body that must work under either writes
    ``Color(getattr(c, 'value', c))``, never ``Color(c)`` and never
    ``c is Color.RED``. See
    ``tests/test_serializer.py::test_enum_arrival_follows_the_node_graph_capability``.

    ``set``/``frozenset`` are deliberately left untouched: a set fails in
    ``general_serializer`` regardless of its contents (not JSON-serializable,
    no registered serializer), so descending into one to flatten enums would
    not make it serializable -- see
    ``tests/test_serializer.py::test_flatten_leaves_sets_untouched``.

    A container with no Enum anywhere inside it comes back as the exact
    object passed in, not a rebuilt plain ``dict``/``tuple``. Rebuilding
    unconditionally would downgrade a ``namedtuple`` to a bare ``tuple`` or
    an ``OrderedDict``/``defaultdict`` to a plain ``dict`` even when nothing
    needed flattening -- see
    ``tests/test_serializer.py::test_flatten_passthrough_preserves_namedtuple_type``
    and ``::test_flatten_passthrough_preserves_dict_subclass_type``.
    """
    if isinstance(value, Enum):
        return _flatten_enums(value.value)
    if isinstance(value, dict):
        out: Dict[Any, Any] = {}
        changed = False
        for k, v in value.items():
            flat_key = _flatten_enums(k)
            flat_val = _flatten_enums(v)
            if flat_key is not k or flat_val is not v:
                changed = True
            if flat_key in out:
                # Two distinct keys (e.g. an Enum member and its bare value,
                # or two members sharing a ``.value``) collapse to one; raise
                # rather than silently drop an entry.
                raise ValueError(f'Enum key flattening collision: multiple keys map to {flat_key!r}')
            out[flat_key] = flat_val
        return out if changed else value
    if isinstance(value, list):
        flat_list = [_flatten_enums(v) for v in value]
        return flat_list if any(fv is not v for fv, v in zip(flat_list, value)) else value
    if isinstance(value, tuple):
        flat_tuple = tuple(_flatten_enums(v) for v in value)
        return flat_tuple if any(fv is not v for fv, v in zip(flat_tuple, value)) else value
    return value


class AiidaSerializationAdapter(SerializationAdapter):
    id: str = 'aiida'
    name: str = 'AiiDA'

    def __init__(self, serializers: Optional[Dict[str, str]] = None, user: Any = None) -> None:
        self.serializers = serializers or all_serializers
        self.user = user

    def serialize(self, value: Any, socket: Any, *, store: bool) -> Any:
        if socket is None:
            return value
        spec = socket._to_spec()
        resolve_tagged_values(value)
        value = _flatten_enums(value)
        return serialize_ports(
            python_data=value,
            port_schema=spec,
            serializers=self.serializers,
            user=self.user,
        )

    def serialize_ports(self, python_data: Any, port_schema: Any, *, store: bool) -> Any:
        resolve_tagged_values(python_data)
        python_data = _flatten_enums(python_data)
        return serialize_ports(
            python_data=python_data,
            port_schema=port_schema,
            serializers=self.serializers,
            user=self.user,
        )
