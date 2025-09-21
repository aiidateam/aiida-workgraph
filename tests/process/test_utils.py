from aiida_workgraph.utils import create_and_pause_process, get_process_summary
from aiida.engine import Process
from aiida import orm
from aiida.common.links import LinkType


class DummyProcess(Process):
    """Name spaced process."""

    _node_class = orm.WorkflowNode


def create_process_node():
    node = orm.CalcJobNode()
    node.set_process_type('aiida.calculations:pythonjob.pythonjob')
    # add inputs
    x = orm.Int(1)
    x.store()
    node.base.links.add_incoming(x, link_type=LinkType.INPUT_CALC, link_label='x')
    node.store()
    # add outputs
    y = orm.Int(1)
    y.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='y')
    y.store()
    return node


def test_create_and_pause_process():
    from aiida.manage import get_manager

    runner = get_manager().get_runner()
    process = create_and_pause_process(runner, DummyProcess, inputs={})
    assert process.node.paused


def test_get_process_summary():
    node = create_process_node()
    summary = get_process_summary(node, data=['inputs', 'outputs'])
    print(summary)
    assert 'Inputs' in summary
    assert 'Outputs' in summary
    assert 'Int' in summary
