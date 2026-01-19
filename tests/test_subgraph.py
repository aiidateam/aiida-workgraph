from aiida_workgraph import WorkGraph, namespace
from typing import Callable
from aiida import orm


def test_inputs_outptus() -> None:
    """Test the inputs and outputs of the WorkGraph."""
    wg = WorkGraph(name='test_inputs_outptus')
    sub_wg = WorkGraph(
        name='sub_wg',
        inputs=namespace(x=int, y=int),
        outputs=namespace(sum=int, diff=int),
    )
    sub_wg_task = wg.add_task(sub_wg, name='sub_wg')
    assert len(sub_wg_task.inputs) == 3
    assert len(sub_wg_task.outputs) == 4
    assert 'x' in sub_wg_task.inputs
    assert 'y' in sub_wg_task.inputs
    assert 'sum' in sub_wg_task.outputs


def test_inputs_outptus_auto_generate(wg_task: WorkGraph) -> None:
    """Test the inputs and outputs of the WorkGraph."""
    wg = WorkGraph(name='test_inputs_outptus')
    # this will generate the group inputs and outputs automatically
    wg_task.expose_inputs()
    wg_task.expose_outputs()
    task1 = wg.add_task(wg_task, name='add1')
    ninput = 0
    for sub_task in wg_task.tasks:
        # remove _wait, but add the namespace
        ninput += len(sub_task.inputs) - 1 + 1
    noutput = 0
    for sub_task in wg_task.tasks:
        noutput += len(sub_task.outputs) - 2 + 1
    assert len(task1.inputs) == 3
    assert len(task1.outputs) == 4
    assert 'sumdiff1.x' in task1.inputs
    assert 'sumdiff1.sum' in task1.outputs


def test_link_subgraph_task(decorated_add: Callable) -> None:
    # create a subgraph
    from node_graph.collection import TaskCollection
    from aiida_workgraph.socket_spec import namespace
    from typing import Any

    x = orm.Int(1).store()

    sub_wg = WorkGraph(
        'build_task_from_workgraph',
        inputs=namespace(x=Any, y=Any),
        outputs=namespace(result=Any),
    )
    sub_wg.add_task(decorated_add, name='add1', x=sub_wg.inputs.x, y=3)
    sub_wg.add_task(
        decorated_add,
        name='add2',
        x=sub_wg.inputs.y,
        y=sub_wg.tasks.add1.outputs.result,
    )
    sub_wg.outputs.result = sub_wg.tasks.add2.outputs.result
    #
    wg = WorkGraph('build_task_from_workgraph')
    add1_task = wg.add_task(decorated_add, name='add1', x=1, y=3)
    wg_task = wg.add_task(sub_wg, name='sub_wg', x=x)
    # the default value of the namespace is None
    assert wg_task.inputs._value == {'x': x}
    assert hasattr(wg.tasks.sub_wg, 'links')
    assert hasattr(wg.tasks.sub_wg, 'tasks')
    assert isinstance(wg.tasks.sub_wg.tasks, TaskCollection)
    assert wg_task.name == 'sub_wg'

    wg.add_task(decorated_add, name='add2', y=3)
    wg.add_link(add1_task.outputs.result, wg_task.inputs.y)
    wg.add_link(wg_task.outputs.result, wg.tasks.add2.inputs.x)
    wg.outputs.sub_wg_result = wg_task.outputs.result
    assert len(wg_task.inputs) == 3
    assert len(wg_task.outputs) == 3
    wg.run()
    assert wg.tasks.add2.outputs.result.value.value == 11
    assert wg.outputs.sub_wg_result.value.value == 8
