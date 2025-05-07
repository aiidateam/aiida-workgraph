from aiida_workgraph import WorkGraph
from typing import Callable
from aiida import orm


def test_inputs_outptus(wg_calcfunction: WorkGraph) -> None:
    """Test the inputs and outputs of the WorkGraph."""
    wg = WorkGraph(name="test_inputs_outptus")
    wg_calcfunction.inputs = {"x": 1, "y": 2}
    wg_calcfunction.outputs.diff = wg_calcfunction.tasks.sumdiff1.outputs.diff
    wg_calcfunction.outputs.sum = wg_calcfunction.tasks.sumdiff2.outputs.sum
    task1 = wg.add_task(wg_calcfunction, name="add1")
    assert len(task1.inputs) == 3
    assert len(task1.outputs) == 4
    assert "x" in task1.inputs
    assert "y" in task1.inputs
    assert "sum" in task1.outputs


def test_inputs_outptus_auto_generate(wg_calcfunction: WorkGraph) -> None:
    """Test the inputs and outputs of the WorkGraph."""
    wg = WorkGraph(name="test_inputs_outptus")
    # this will generate the group inputs and outputs automatically
    task1 = wg.add_task(wg_calcfunction, name="add1")
    ninput = 0
    for sub_task in wg_calcfunction.tasks:
        # remove _wait, but add the namespace
        ninput += len(sub_task.inputs) - 1 + 1
    noutput = 0
    for sub_task in wg_calcfunction.tasks:
        noutput += len(sub_task.outputs) - 2 + 1
    assert len(task1.inputs) == len(wg_calcfunction.tasks) + 1
    assert len(task1.outputs) == len(wg_calcfunction.tasks) + 2
    assert "sumdiff1.x" in task1.inputs
    assert "sumdiff1.sum" in task1.outputs


def test_build_task_from_workgraph(decorated_add: Callable) -> None:
    # create a sub workgraph
    from aiida_workgraph.collection import TaskCollection

    x = orm.Int(1).store()

    sub_wg = WorkGraph("build_task_from_workgraph")
    sub_wg.add_task(decorated_add, name="add1", x=x, y=3)
    sub_wg.add_task(decorated_add, name="add2", x=2, y=sub_wg.tasks.add1.outputs.result)
    #
    wg = WorkGraph("build_task_from_workgraph")
    add1_task = wg.add_task(decorated_add, name="add1", x=1, y=3)
    wg_task = wg.add_task(sub_wg, name="sub_wg")
    # the default value of the namespace is None
    assert wg_task.inputs["add1"]._value == {"x": x, "y": 3}
    assert hasattr(wg.tasks.sub_wg, "workgraph")
    assert hasattr(wg.tasks.sub_wg, "links")
    assert hasattr(wg.tasks.sub_wg, "tasks")
    assert isinstance(wg.tasks.sub_wg.tasks, TaskCollection)
    assert wg.tasks.sub_wg.tasks.add1.graph.name == "build_task_from_workgraph"

    wg.add_task(decorated_add, name="add2", y=3)
    wg.add_link(add1_task.outputs.result, wg_task.inputs["add1.x"])
    wg.add_link(wg_task.outputs["add2.result"], wg.tasks.add2.inputs.x)
    assert len(wg_task.inputs) == 3
    assert len(wg_task.outputs) == 4
    wg.run()
    assert wg.tasks.add2.outputs.result.value.value == 12
