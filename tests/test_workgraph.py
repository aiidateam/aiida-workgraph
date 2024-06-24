from aiida_workgraph import WorkGraph, build_task
from aiida import load_profile, orm
import time
import pytest
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

load_profile()


def test_to_dict(wg_calcjob):
    """Export NodeGraph to dict."""
    wg = wg_calcjob
    wgdata = wg.to_dict()
    assert len(wgdata["tasks"]) == len(wg.tasks)
    assert len(wgdata["links"]) == len(wg.links)


def test_from_dict(wg_calcjob):
    """Export NodeGraph to dict."""
    wg = wg_calcjob
    wgdata = wg.to_dict()
    wg1 = WorkGraph.from_dict(wgdata)
    assert len(wg.tasks) == len(wg1.tasks)
    assert len(wg.links) == len(wg1.links)


def test_new_node(wg_calcjob):
    """Add new task."""
    wg = wg_calcjob
    n = len(wg.tasks)
    wg.tasks.new(ArithmeticAddCalculation)
    assert len(wg.tasks) == n + 1


def test_save_load(wg_calcjob):
    """Save the workgraph"""
    wg = wg_calcjob
    wg.name = "test_save_load"
    wg.save()
    assert wg.process.process_state.value.upper() == "CREATED"
    assert wg.process.process_label == "WorkGraph<test_save_load>"
    wg2 = WorkGraph.load(wg.process.pk)
    assert len(wg.tasks) == len(wg2.tasks)


# skip this test
@pytest.mark.skip(reason="PAUSED state is wrong for the moment.")
def test_pause(wg_engine):
    wg = wg_engine
    wg.name = "test_pause"
    wg.submit()
    time.sleep(5)
    wg.pause()
    wg.update()
    assert wg.process.process_state.value.upper() == "PAUSED"


def test_reset_message(wg_calcjob):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    from aiida.cmdline.utils.common import get_workchain_report

    wg = wg_calcjob
    wg.submit()
    wg = WorkGraph.load(wg.process.pk)
    wg.tasks["add2"].set({"y": orm.Int(10).store()})
    wg.save()
    wg.wait()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Task add2 action: RESET." in report


def test_restart(wg_calcjob):
    """Restart from a finished workgraph.
    Load the workgraph, modify the task, and restart the workgraph.
    Only the modified node and its child tasks will be rerun."""
    wg = wg_calcjob
    wg.name = "test_restart_0"
    wg.submit(wait=True)
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.name = "test_restart_1"
    wg1.tasks["add2"].set({"y": orm.Int(10).store()})
    wg1.submit(wait=True)
    wg1.update()
    assert wg1.tasks["add3"].node.outputs.sum == 13
    assert wg1.tasks["add1"].node.pk == wg.tasks["add1"].pk
    assert wg1.tasks["add2"].node.pk != wg.tasks["add2"].pk


def test_extend_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph("test_graph_build")
    add1 = wg.tasks.new("AiiDAAdd", "add1", x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    # extend workgraph
    wg.extend(add_multiply_wg, prefix="group_")
    wg.links.new(add1.outputs[0], wg.tasks["group_add1"].inputs["x"])
    wg.submit(wait=True)
    assert wg.tasks["group_multiply1"].node.outputs.result == 45


def test_node_from_workgraph(decorated_add_multiply_group):
    wg = WorkGraph("test_node_from_workgraph")
    add1 = wg.tasks.new("AiiDAAdd", "add1", x=2, y=3)
    add2 = wg.tasks.new("AiiDAAdd", "add2", y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    AddMultiplyTask = build_task(add_multiply_wg)
    assert "add1.x" in AddMultiplyTask().inputs.keys()
    # add the workgraph as a task
    add_multiply1 = wg.tasks.new(AddMultiplyTask, "add_multiply1")
    wg.links.new(add1.outputs[0], add_multiply1.inputs["add1.x"])
    wg.links.new(add_multiply1.outputs["multiply1.result"], add2.inputs["x"])
    # wg.submit(wait=True)
    wg.run()
    assert wg.tasks["add2"].node.outputs.sum == 48


def test_pause_task(wg_calcjob):
    wg = wg_calcjob
    wg.name = "test_pause_task"
    wg.submit()
    # wg.run()
    wg.pause_tasks(["add2"])
    time.sleep(20)
    wg.update()
    assert wg.tasks["add2"].node.process_state.value.upper() == "CREATED"
    assert wg.tasks["add2"].node.process_status == "Paused through WorkGraph"
    wg.play_tasks(["add2"])
    wg.wait()
    assert wg.tasks["add2"].outputs["sum"].value == 9
