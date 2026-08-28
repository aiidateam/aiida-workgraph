from collections import OrderedDict, defaultdict, namedtuple
from enum import Enum, IntEnum

from aiida_workgraph import WorkGraph, task
import pytest

from aiida_workgraph.serialization import AiidaSerializationAdapter, _flatten_enums


@task.graph()
def sub_workflow(func):
    func()


def test_func_as_input(capsys):
    from aiida_workgraph.executors.test import add

    wg = WorkGraph('test_func_as_input')
    wg.add_task(sub_workflow, func=add, name='sub_workflow')
    with pytest.raises(Exception, match='Cannot serialize the provided object'):
        wg.save()


class Color(Enum):
    """Plain enum whose value is a bare string."""

    RED = 'red'
    BLUE = 'blue'


class Priority(IntEnum):
    """Enum whose value is a bare int."""

    LOW = 1
    HIGH = 9


class Flavour(str, Enum):
    """Enum member that is itself a ``str`` subclass."""

    SWEET = 'sweet'
    SOUR = 'sour'


def test_flatten_bare_enum_member():
    """A lone member collapses to its ``.value``, not a mangled ``Color.RED``."""
    assert _flatten_enums(Color.RED) == 'red'
    assert not isinstance(_flatten_enums(Color.RED), Enum)


def test_flatten_intenum_and_str_enum_yield_bare_value():
    """The bare value's own type is preserved (int stays int, str stays str)."""
    assert _flatten_enums(Priority.HIGH) == 9
    assert type(_flatten_enums(Priority.HIGH)) is int
    assert _flatten_enums(Flavour.SWEET) == 'sweet'
    assert type(_flatten_enums(Flavour.SWEET)) is str


def test_flatten_dict_keys_and_values():
    """Enums appearing as dict keys *and* values both flatten."""
    assert _flatten_enums({Color.RED: Color.BLUE, 'n': Priority.LOW}) == {'red': 'blue', 'n': 1}


def test_flatten_list_and_tuple_preserve_container_type():
    """Sequences recurse element-wise and keep list-vs-tuple identity."""
    assert _flatten_enums([Color.RED, 1, 'x']) == ['red', 1, 'x']
    out = _flatten_enums((Color.RED, Priority.HIGH))
    assert out == ('red', 9)
    assert isinstance(out, tuple)


def test_flatten_mixed_nesting():
    """Enums buried in a dict/list/tuple mixture are all reached."""
    payload = {
        'a': [Color.RED, {'b': Priority.LOW}],
        Color.BLUE: ('x', Flavour.SOUR),
    }
    assert _flatten_enums(payload) == {
        'a': ['red', {'b': 1}],
        'blue': ('x', 'sour'),
    }


def test_flatten_passthrough_of_enum_free_payload():
    """An enum-free payload comes back with identical values."""
    payload = {'n': 3, 'items': [1, 'two', (3.0, None)], 'flag': True}
    assert _flatten_enums(payload) == payload


def test_flatten_passthrough_preserves_namedtuple_type():
    """An enum-free namedtuple comes back as the exact same object, not
    downgraded to a plain tuple."""
    Point = namedtuple('Point', ['x', 'y'])
    p = Point(1, 2)
    out = _flatten_enums(p)
    assert out is p
    assert isinstance(out, Point)


def test_flatten_namedtuple_with_enum_still_flattens():
    """A namedtuple carrying an Enum member still flattens its contents."""
    Point = namedtuple('Point', ['x', 'y'])
    p = Point(Color.RED, 2)
    out = _flatten_enums(p)
    assert out == ('red', 2)


def test_flatten_passthrough_preserves_dict_subclass_type():
    """An enum-free OrderedDict/defaultdict comes back as the exact same
    object, not downgraded to a plain dict."""
    od = OrderedDict([('a', 1), ('b', 2)])
    assert _flatten_enums(od) is od

    dd: defaultdict = defaultdict(int, {'a': 1})
    assert _flatten_enums(dd) is dd


def test_flatten_dict_subclass_with_enum_still_flattens():
    """An OrderedDict/defaultdict carrying an Enum still flattens its
    values (into a plain dict -- rebuilding the original subclass is not
    attempted for the changed case)."""
    od = OrderedDict([('a', Color.RED), ('b', 2)])
    assert _flatten_enums(od) == {'a': 'red', 'b': 2}

    dd: defaultdict = defaultdict(int, {'a': Priority.LOW})
    assert _flatten_enums(dd) == {'a': 1}


def test_flatten_wrapt_proxied_enum_member():
    """A wrapt ``ObjectProxy`` (as TaggedValue wraps members) still flattens."""
    wrapt = pytest.importorskip('wrapt')
    proxied = wrapt.ObjectProxy(Color.RED)
    assert _flatten_enums(proxied) == 'red'
    assert _flatten_enums({proxied: proxied}) == {'red': 'red'}


def test_flatten_dict_key_collision_raises():
    """Two distinct keys collapsing to the same flattened key raises rather than
    silently dropping an entry. Here ``Color.RED`` (a plain Enum, so distinct in
    identity and hash from the bare string) and ``'red'`` are separate keys that
    both flatten to ``'red'``."""
    payload = {Color.RED: 1, 'red': 2}
    assert len(payload) == 2  # distinct keys before flattening
    with pytest.raises(ValueError, match='collision'):
        _flatten_enums(payload)


def test_flatten_leaves_sets_untouched():
    """``set``/``frozenset`` are not descended into: they fail in
    ``general_serializer`` regardless of contents, so flattening enums inside
    them would not make them serializable. Pin that they pass through unchanged
    (still holding Enum members)."""
    s = {Color.RED, Color.BLUE}
    out = _flatten_enums(s)
    assert out is s
    assert Color.RED in out
    fs = frozenset({Priority.LOW})
    assert _flatten_enums(fs) is fs


def test_serialize_ports_accepts_enum_value(aiida_profile):
    """End-to-end through the adapter: an enum-valued namespace entry serializes
    without raising, landing as the flattened bare value (needs a profile so the
    downstream ``general_serializer`` can build the node)."""
    from aiida import orm
    from node_graph.socket_spec import SocketMeta, SocketSpec

    spec = SocketSpec(identifier='node_graph.namespace', meta=SocketMeta(dynamic=True))
    out = AiidaSerializationAdapter().serialize_ports({'c': Color.RED, 'n': 3}, spec, store=False)
    assert isinstance(out['c'], orm.Str)
    assert out['c'].value == 'red'


@task(outputs=['type_name', 'is_member', 'equals_member', 'equals_value'])
def observe_enum(c: Color) -> dict:
    """Report what a function task's body receives for an ``Enum``-typed
    input, using only the annotation's contract (``==``), never a
    reconstruction of the member."""
    return {
        'type_name': type(c).__name__,
        'is_member': isinstance(c, Color),
        'equals_member': c == Color.RED,
        'equals_value': c == 'red',
    }


@task()
def echo_str(s: str) -> str:
    return s


@task.graph()
def observe_enum_in_graph(c: Color) -> str:
    """Use the member a ``@task.graph`` body receives, as the annotation
    declares it — no reconstruction."""
    return echo_str(s=c.name)


def _node_graph_rebuilds_enums() -> bool:
    """Return whether the installed ``node_graph`` reconstructs Enum sockets.

    ``coerce_inputs_from_spec`` rebuilds a socket's value only from the
    ``structured_type`` descriptor its spec records, so a descriptor for an
    ``Enum`` type is exactly the capability.
    """
    from node_graph.utils.struct_utils import structured_type_info

    return structured_type_info(Color) is not None


def test_enum_arrival_follows_the_node_graph_capability(aiida_profile):
    """A task body never reconstructs an ``Enum`` member by hand; the
    annotation is the contract, and what arrives behind it is the installed
    ``node_graph``'s call, not this package's. With a ``node_graph`` that
    rebuilds ``Enum`` sockets (``_node_graph_rebuilds_enums``), a
    ``Color``-annotated input arrives as the member: ``isinstance`` holds,
    ``==`` against the member holds, ``type_name`` is ``'Color'``. Without
    that capability it arrives as the flattened bare value: ``isinstance``
    fails, ``==`` against the value holds, ``type_name`` is ``'str'``."""
    wg = WorkGraph('enum_arrival')
    t = wg.add_task(observe_enum, name='observe', c=Color.RED)
    wg.run()
    assert wg.state == 'FINISHED'
    rebuilds = _node_graph_rebuilds_enums()
    assert t.outputs['is_member'].value.value is rebuilds
    assert t.outputs['type_name'].value.value == ('Color' if rebuilds else 'str')
    if rebuilds:
        assert t.outputs['equals_member'].value.value is True
    else:
        assert t.outputs['equals_value'].value.value is True


def test_enum_input_arrives_as_the_member_in_a_graph_body(aiida_profile):
    """Only under a ``node_graph`` that rebuilds ``Enum`` sockets
    (``_node_graph_rebuilds_enums``) can a ``@task.graph`` body use the
    member directly, as its annotation declares (``c.name``, no
    reconstruction)."""
    if not _node_graph_rebuilds_enums():
        pytest.skip(
            'node_graph does not reconstruct Enum sockets '
            '(structured_type_info(Color) is None); '
            'observe_enum_in_graph requires the member to use c.name.'
        )
    wg = WorkGraph('enum_arrival_in_graph')
    t = wg.add_task(observe_enum_in_graph, name='observe_graph', c=Color.RED)
    wg.run()
    assert wg.state == 'FINISHED'
    assert t.outputs['result'].value.value == 'RED'
