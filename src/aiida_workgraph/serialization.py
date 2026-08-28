from __future__ import annotations

from typing import Any, Dict, Optional

from aiida_pythonjob.data.serializer import all_serializers
from aiida_pythonjob.utils import serialize_ports
from node_graph.serializer import SerializationAdapter
from node_graph.utils import resolve_tagged_values

# Socket identifiers whose declared type is a Python primitive. When the
# persisted value round-trips through AiiDA storage it becomes the matching
# ``orm.BaseType`` node (``orm.Float``, ``orm.Int``, ``orm.Str``,
# ``orm.Bool``); ``deserialize`` must strip that wrapper before the value
# reaches a user-written ``@task.graph`` body whose signature declared a
# primitive type.
_PRIMITIVE_SOCKET_IDENTIFIERS = frozenset(
    {
        'workgraph.float',
        'workgraph.int',
        'workgraph.string',
        'workgraph.bool',
    }
)


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
        return serialize_ports(
            python_data=value,
            port_schema=spec,
            serializers=self.serializers,
            user=self.user,
        )

    def deserialize(self, value: Any, socket: Any) -> Any:
        """Unwrap ``orm.BaseType`` to its Python value for primitive sockets.

        Symmetric counterpart to ``serialize``: where the write path
        auto-serialises primitive inputs into ``orm.Float`` / ``orm.Int`` /
        ``orm.Str`` / ``orm.Bool`` for provenance, the read path (invoked
        just before a ``@task.graph`` body is called) hands the body the
        primitive that its signature actually declared. Non-primitive sockets
        (``workgraph.any``, AiiDA-typed sockets, etc.) pass through.

        Dataclass-typed sockets get the same unwrap applied to their
        fields. Without this the ``cls(**value)`` reconstruction in
        ``node_graph.utils.struct_utils.coerce_structured_value`` re-packs
        node-promoted scalars into the dataclass, so a field declared
        ``int`` lands as ``orm.Int`` and downstream stdlib calls
        (``range(self.ntyp)``, etc.) fail with ``TypeError`` because
        ``orm.Int`` doesn't implement ``__index__``.

        A ``TaggedValue`` -- node-graph's proxy marking a value that still
        draws a link back to its owning graph-input socket -- keeps its tag:
        the unwrap runs on the value it wraps, and the result is rewrapped in
        a fresh ``TaggedValue`` pointing at the same socket. Returning the
        bare unwrapped value instead makes a sub-task bound from it store a
        new, unlinked copy of the graph input rather than draw provenance
        from it.
        """
        from dataclasses import fields, is_dataclass, replace

        from aiida import orm
        from node_graph.socket import TaggedValue

        if isinstance(value, TaggedValue):
            tag = value._socket
            unwrapped = self.deserialize(value.__wrapped__, socket)
            if unwrapped is value.__wrapped__:
                return value
            return TaggedValue(unwrapped, socket=tag)

        if isinstance(value, orm.BaseType):
            identifier = getattr(socket, '_identifier', None)
            if identifier in _PRIMITIVE_SOCKET_IDENTIFIERS:
                return value.value
            return value

        # Same unwrap for container-typed sockets: a plain dict/list built in
        # a parent graph body is promoted to ``orm.Dict`` / ``orm.List`` for
        # provenance, but a body whose signature declares ``dict`` / ``list``
        # expects the plain container (e.g. ``get_protocol_inputs`` calls
        # ``overrides.copy()``, which ``orm.Dict`` doesn't implement).
        identifier = getattr(socket, '_identifier', None)
        if isinstance(value, orm.Dict) and identifier == 'workgraph.dict':
            return value.get_dict()
        if isinstance(value, orm.List) and identifier == 'workgraph.list':
            return value.get_list()

        if is_dataclass(value) and not isinstance(value, type):
            # Each field is tagged independently of the dataclass instance as
            # a whole (``tag_socket_value`` walks a structured socket down to
            # its leaves), so the tag to preserve lives on the field value,
            # not on ``value`` itself.
            field_updates = {}
            for f in fields(value):
                field_value = getattr(value, f.name)
                tag = None
                raw = field_value
                if isinstance(field_value, TaggedValue):
                    tag = field_value._socket
                    raw = field_value.__wrapped__
                if isinstance(raw, orm.BaseType):
                    plain = raw.value
                    field_updates[f.name] = TaggedValue(plain, socket=tag) if tag is not None else plain
            if field_updates:
                return replace(value, **field_updates)

        return value

    def serialize_ports(self, python_data: Any, port_schema: Any, *, store: bool) -> Any:
        resolve_tagged_values(python_data)
        return serialize_ports(
            python_data=python_data,
            port_schema=port_schema,
            serializers=self.serializers,
            user=self.user,
        )
