import enum
import json

import pytest
from aiida import orm
from aiida.common.exceptions import ValidationError

from aiida_workgraph.utils import _ensure_json_safe, _ensure_json_safe_key


class _Color(enum.Enum):
    """Plain Enum: NOT a str/int subclass -> not JSON-serializable."""

    RED = 1
    BLUE = 2


class _SpinChannel(str, enum.Enum):
    """str-Enum: IS a str instance, so already JSON-serializable."""

    UP = 'up'
    DOWN = 'down'


class _Count(enum.IntEnum):
    """IntEnum: IS an int instance, so already JSON-serializable."""

    ONE = 1


def test_ensure_json_safe_plain_enum_unwraps_to_value():
    """A plain Enum (the real bug class) is unwrapped to its ``.value``."""
    assert _ensure_json_safe(_Color.RED) == 1
    assert _ensure_json_safe(_Color.BLUE) == 2
    # nested exactly as an enum default would sit inside wgdata
    wgdata = {'tasks': {'t': {'inputs': {'color': {'value': _Color.RED}}}}}
    assert _ensure_json_safe(wgdata) == {'tasks': {'t': {'inputs': {'color': {'value': 1}}}}}


def test_ensure_json_safe_str_enum_is_noop():
    """A str-Enum is already a str, so the helper leaves it untouched.

    This documents that the helper is a no-op for str-Enum: it must NOT be used
    as a regression fixture for the bug, since it passes with or without the
    helper (see PR discussion).
    """
    assert _ensure_json_safe(_SpinChannel.UP) is _SpinChannel.UP
    assert json.dumps(_ensure_json_safe(_SpinChannel.UP)) == '"up"'


def test_ensure_json_safe_int_enum_serializes():
    """An IntEnum is already an int; unchanged and JSON-serializable."""
    assert _ensure_json_safe(_Count.ONE) is _Count.ONE
    assert json.dumps(_ensure_json_safe(_Count.ONE)) == '1'


def test_ensure_json_safe_dict_enum_key_is_coerced():
    """Enum dict keys are coerced (regression: keys were previously ignored)."""
    result = _ensure_json_safe({_Color.RED: 'x', _Color.BLUE: 'y'})
    assert result == {1: 'x', 2: 'y'}
    # the whole thing must now be JSON-serializable (keys included)
    assert json.dumps(result) == '{"1": "x", "2": "y"}'


def test_ensure_json_safe_non_enum_key_stringified():
    """A non-enum, non-primitive key is stringified so json.dumps succeeds."""

    class _Key:
        def __str__(self):
            return 'k'

    assert _ensure_json_safe_key(_Key()) == 'k'
    assert _ensure_json_safe({_Key(): 1}) == {'k': 1}


def test_ensure_json_safe_preserves_clean_value_coercible():
    """Values that ``clean_value`` coerces itself pass through untouched.

    The helper must not pre-empt storage's own coercion: a set becomes a list
    at store time, a numpy scalar a Python scalar, a ``BaseType`` its value —
    none of them may be stringified by the helper.
    """
    import numpy as np

    assert _ensure_json_safe({1, 2, 3}) == {1, 2, 3}
    assert _ensure_json_safe(frozenset({1, 2})) == frozenset({1, 2})
    numpy_int = np.int64(7)
    assert _ensure_json_safe(numpy_int) is numpy_int
    aiida_int = orm.Int(3)
    assert _ensure_json_safe(aiida_int) is aiida_int

    node = orm.WorkflowNode()
    node.base.attributes.set('data', _ensure_json_safe({'tags': {1, 2, 3}, 'n': np.int64(7)}))
    node.store()
    stored = orm.load_node(node.pk).base.attributes.get('data')
    assert sorted(stored['tags']) == [1, 2, 3]
    assert stored['n'] == 7


def test_ensure_json_safe_non_dict_mapping_is_coerced():
    """A non-``dict`` ``Mapping`` takes the mapping branch, keys included.

    ``clean_value`` does not inspect mapping keys, so a bad key inside e.g. a
    ``MappingProxyType`` would otherwise pass the helper and only fail in the
    database driver at store time.
    """
    from types import MappingProxyType

    result = _ensure_json_safe(MappingProxyType({(1, 2): 'v', _Color.RED: 'w'}))
    assert result == {'(1, 2)': 'v', 1: 'w'}

    node = orm.WorkflowNode()
    node.base.attributes.set('data', result)
    node.store()
    assert orm.load_node(node.pk).base.attributes.get('data') == {'(1, 2)': 'v', '1': 'w'}


def test_ensure_json_safe_iterator_is_materialized():
    """A one-shot iterator is materialized to a list, not passed to
    ``clean_value`` (which would exhaust it as a side effect of validation and
    leave an empty value to be stored)."""
    assert _ensure_json_safe(iter([1, 2])) == [1, 2]
    assert _ensure_json_safe(v for v in (_Color.RED, 3)) == [1, 3]


def test_ensure_json_safe_fallback_does_not_unwrap_arbitrary_value_attr():
    """The fallback is narrowed to Enum: a plain object exposing ``.value`` is
    NOT silently unwrapped (it is stringified instead)."""

    class _Wrapper:
        value = 42

    result = _ensure_json_safe(_Wrapper())
    assert result != 42
    assert isinstance(result, str)
    json.dumps(result)  # must not raise


def test_ensure_json_safe_output_is_storable():
    """Whatever the helper returns must always survive ``clean_value``."""
    from aiida.orm.implementation.utils import clean_value

    payload = {
        'a_plain_enum': _Color.RED,
        'a_str_enum': _SpinChannel.DOWN,
        'a_list': [_Color.BLUE, (1, 2), {_Color.RED: 'nested'}],
        'a_set': {1, 2},
        _Color.RED: 'enum-key',
    }
    clean_value(_ensure_json_safe(payload))  # must not raise


def test_enum_attribute_store_negative_control():
    """End-to-end at the real storage gate (``clean_value`` + JSONB).

    Negative control: a plain-Enum-bearing dict set as a node attribute raises
    ``ValidationError`` when the node is stored (``clean_value`` runs at store
    time); wrapping it in ``_ensure_json_safe`` first lets the store succeed and
    the value round-trips through the database.
    """
    wgdata = {'tasks': {'t': {'inputs': {'color': {'value': _Color.RED}}}}}

    raw_node = orm.WorkflowNode()
    raw_node.base.attributes.set('workgraph_data', wgdata)
    with pytest.raises(ValidationError):
        raw_node.store()

    safe_node = orm.WorkflowNode()
    safe_node.base.attributes.set('workgraph_data', _ensure_json_safe(wgdata))
    safe_node.store()
    stored = orm.load_node(safe_node.pk).base.attributes.get('workgraph_data')
    assert stored == {'tasks': {'t': {'inputs': {'color': {'value': 1}}}}}


def test_save_coerces_graph_level_error_handler_kwargs():
    """`WorkGraph.save()` must coerce graph-level error-handler ``kwargs``.

    Error-handler ``kwargs`` are copied verbatim into the
    ``workgraph_error_handlers`` attribute, so a plain-Enum value there reaches
    the storage gate.  This exercises the wiring of ``_ensure_json_safe`` inside
    ``save_workgraph_data``, not just the helper in isolation.
    """
    from node_graph.error_handler import normalize_error_handlers

    from aiida_workgraph import WorkGraph, task

    @task()
    def add(x: int = 1, y: int = 2):
        return x + y

    def handle(task, **kwargs):  # never runs, only stored
        return 'retrying'

    wg = WorkGraph(
        'graph_level_handler',
        error_handlers=normalize_error_handlers(
            {'h': {'executor': handle, 'exit_codes': [1], 'kwargs': {'color': _Color.RED}}}
        ),
    )
    wg.add_task(add, name='add1')
    wg.save()

    stored = orm.load_node(wg.process.pk).base.attributes.get('workgraph_error_handlers')
    assert stored['h']['kwargs']['color'] == 1


def test_save_coerces_task_level_error_handler_kwargs():
    """`WorkGraph.save()` must coerce task-level error-handler ``kwargs``.

    A task-level handler is stored inside the ``workgraph_data`` attribute under
    that task's spec; a plain-Enum value in its ``kwargs`` must be unwrapped
    there too.
    """
    from aiida_workgraph import WorkGraph, task

    @task()
    def add(x: int = 1, y: int = 2):
        return x + y

    def handle(task, **kwargs):  # never runs, only stored
        return 'retrying'

    wg = WorkGraph('task_level_handler')
    add1 = wg.add_task(add, name='add1')
    add1.add_error_handler({'h': {'executor': handle, 'exit_codes': [1], 'kwargs': {'color': _Color.RED}}})
    wg.save()

    stored = orm.load_node(wg.process.pk).base.attributes.get('workgraph_data')
    handlers = stored['tasks']['add1']['spec']['attached_error_handlers']
    assert handlers['h']['kwargs']['color'] == 1


def test_get_or_create_code(fixture_localhost):
    from aiida_workgraph.utils import get_or_create_code
    from aiida.orm import Code

    # create a new code
    code1 = get_or_create_code(
        computer='localhost',
        code_label='test_code',
        code_path='/bin/bash',
        prepend_text='echo "Hello, World!"',
    )
    assert isinstance(code1, Code)
    # use already created code
    code2 = get_or_create_code(
        computer='localhost',
        code_label='test_code',
        code_path='/bin/bash',
        prepend_text='echo "Hello, World!"',
    )
    assert code1.uuid == code2.uuid


def test_get_parent_workgraphs():
    from aiida.common.links import LinkType
    from aiida_workgraph.utils import get_parent_workgraphs

    wn1 = orm.WorkflowNode()
    wn2 = orm.WorkflowNode()
    wn3 = orm.WorkflowNode()
    wn3.base.links.add_incoming(wn2, link_type=LinkType.CALL_WORK, link_label='link')
    wn2.base.links.add_incoming(wn1, link_type=LinkType.CALL_WORK, link_label='link')
    wn1.store()
    wn2.store()
    wn3.store()

    parent_workgraphs = get_parent_workgraphs(wn3.pk)
    assert len(parent_workgraphs) == 3


def test_generate_provenance_graph():
    from IPython.display import IFrame
    from aiida_workgraph.utils import generate_provenance_graph
    import os

    wn1 = orm.WorkflowNode()
    wn1.store()

    graph = generate_provenance_graph(wn1.pk)
    assert isinstance(graph, IFrame)
    # check file html/node_graph_{pk}.html is created
    assert os.path.isfile(f'html/node_graph_{wn1.pk}.html')


def test_workgraph_generate_provenance_graph():
    from IPython.display import IFrame
    from aiida_workgraph import WorkGraph
    import os

    wg = WorkGraph()
    wg.save()

    graph = wg.generate_provenance_graph()
    assert isinstance(graph, IFrame)
    assert os.path.isfile(f'html/node_graph_{wg.pk}.html')
