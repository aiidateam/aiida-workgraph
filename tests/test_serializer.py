from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from aiida_workgraph import WorkGraph, task
from aiida_workgraph.serialization import AiidaSerializationAdapter
import pytest


@task.graph()
def sub_workflow(func):
    func()


def test_func_as_input(capsys):
    from aiida_workgraph.executors.test import add

    wg = WorkGraph('test_func_as_input')
    wg.add_task(sub_workflow, func=add, name='sub_workflow')
    with pytest.raises(Exception, match='Cannot serialize the provided object'):
        wg.save()


# --- deserialize must not break provenance while unwrapping ------------------
#
# A ``@task.graph`` body whose signature declares a primitive type receives
# the unwrapped ``str``/``int``/``dict``/``list``/dataclass, not the
# ``orm.BaseType``/``orm.Dict``/``orm.List`` node the write path promoted it
# to. If a sub-task inside the body is then bound from that value, the value
# must still carry the ``TaggedValue`` tag that draws a link back to the
# graph input's own node -- otherwise the sub-task stores a brand-new,
# unlinked copy (an orphan input with zero incoming links) instead of
# referencing the graph input.


@task()
def echo_task(label, count):
    return f'{label}-{count}'


@task.graph()
def bare_graph(label, count):
    return echo_task(label=label, count=count)


@task.graph()
def annotated_graph(label: str, count: int):
    return echo_task(label=label, count=count)


@task.graph()
def annotated_any_graph(label: Any, count: Any):
    return echo_task(label=label, count=count)


def _label_link_uuids(node):
    """Return (graph-input node, sub-task input node) for the ``label`` link."""
    outer = inner = None
    for link in node.process.base.links.get_incoming().all():
        if link.link_label.endswith('label'):
            outer = link.node
    for child in node.process.called:
        for link in child.base.links.get_incoming().all():
            if link.link_label.endswith('label'):
                inner = link.node
    return outer, inner


@pytest.mark.parametrize(
    'entry, name', [(bare_graph, 'bare'), (annotated_graph, 'annotated'), (annotated_any_graph, 'annotated_any')]
)
def test_primitive_graph_input_link_survives_deserialize(entry, name):
    """A primitive graph input still links to the sub-task it feeds.

    ``annotated`` is the case #799 broke: ``label: str`` makes ``deserialize``
    unwrap the socket's ``orm.Str`` to a bare ``str`` before the body runs, and
    without tag preservation that bare ``str`` becomes a fresh, unlinked node
    when ``echo_task`` is bound from it.
    """
    wg = WorkGraph(f'link_{name}')
    node = wg.add_task(entry, name='g', label='silicon', count=2)
    wg.run()
    outer, inner = _label_link_uuids(node)
    assert outer is not None and inner is not None
    assert outer.uuid == inner.uuid, (
        f'sub-task input {inner.uuid} is an orphan copy of graph input {outer.uuid}, '
        f'not a link to it (incoming links: {len(inner.base.links.get_incoming().all())})'
    )


@task()
def echo_container_task(d, lst):
    return {'d': d, 'l': lst}


@task.graph()
def container_graph(d: dict, lst: list):
    return echo_container_task(d=d, lst=lst)


def test_dict_and_list_graph_input_link_survives_deserialize():
    """A dict/list graph input, unwrapped from ``orm.Dict``/``orm.List``, still links."""
    wg = WorkGraph('link_container')
    node = wg.add_task(container_graph, name='g', d={'a': 1}, lst=[1, 2])
    wg.run()
    outer = inner = None
    for link in node.process.base.links.get_incoming().all():
        if link.link_label.endswith('d'):
            outer = link.node
    for child in node.process.called:
        for link in child.base.links.get_incoming().all():
            if link.link_label.endswith('d'):
                inner = link.node
    assert outer is not None and inner is not None
    assert outer.uuid == inner.uuid, 'dict graph input was copied instead of linked'


@dataclass
class _Named:
    label: str
    count: int


@task.graph()
def dataclass_graph(data: _Named):
    return echo_task(label=data.label, count=data.count)


def test_dataclass_field_link_survives_deserialize():
    """A dataclass-typed graph input's fields are unwrapped without losing their links.

    Each field is tagged independently of the dataclass instance (a
    structured socket is tagged leaf by leaf), so the fix must preserve the
    tag on the field value it unwraps, not on the dataclass instance itself.
    """
    wg = WorkGraph('link_dataclass')
    node = wg.add_task(dataclass_graph, name='g', data=_Named(label='silicon', count=2))
    wg.run()
    outer, inner = _label_link_uuids(node)
    assert outer is not None and inner is not None
    assert outer.uuid == inner.uuid, 'dataclass field was copied instead of linked'


def _pre_fix_deserialize(self, value, socket):
    """The #799 body before this fix: unwraps via attribute access with no
    awareness of ``TaggedValue``, so the tag that draws a provenance link is
    silently dropped. Used as a negative control below."""
    from dataclasses import fields, is_dataclass, replace

    from aiida import orm

    if isinstance(value, orm.BaseType):
        identifier = getattr(socket, '_identifier', None)
        if identifier in {'workgraph.float', 'workgraph.int', 'workgraph.string', 'workgraph.bool'}:
            return value.value
        return value

    identifier = getattr(socket, '_identifier', None)
    if isinstance(value, orm.Dict) and identifier == 'workgraph.dict':
        return value.get_dict()
    if isinstance(value, orm.List) and identifier == 'workgraph.list':
        return value.get_list()

    if is_dataclass(value) and not isinstance(value, type):
        field_updates = {
            f.name: getattr(value, f.name).value
            for f in fields(value)
            if isinstance(getattr(value, f.name), orm.BaseType)
        }
        if field_updates:
            return replace(value, **field_updates)

    return value


def test_negative_control_pre_fix_deserialize_loses_the_link(monkeypatch):
    """Discriminating control: the pre-fix unwrap reproduces #799's link loss.

    If this failed to reproduce the loss, the provenance assertions above
    would prove nothing -- they could be passing for a reason unrelated to
    the fix.
    """
    monkeypatch.setattr(AiidaSerializationAdapter, 'deserialize', _pre_fix_deserialize)
    wg = WorkGraph('link_negative_control')
    node = wg.add_task(annotated_graph, name='g', label='silicon', count=2)
    wg.run()
    outer, inner = _label_link_uuids(node)
    assert outer is not None and inner is not None
    assert outer.uuid != inner.uuid, 'pre-fix deserialize unexpectedly preserved the link'
    assert len(inner.base.links.get_incoming().all()) == 0, 'expected an orphan copy'
