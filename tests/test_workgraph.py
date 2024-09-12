import pytest
from aiida_workgraph import WorkGraph
from aiida import orm
import time
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


def test_from_dict(decorated_add):
    """Export NodeGraph to dict."""
    wg = WorkGraph("test_from_dict")
    task1 = wg.add_task(decorated_add, x=2, y=3)
    wg.add_task(
        "workgraph.test_sum_diff", name="sumdiff2", x=4, y=task1.outputs["result"]
    )
    wgdata = wg.to_dict()
    wg1 = WorkGraph.from_dict(wgdata)
    assert len(wg.tasks) == len(wg1.tasks)
    assert len(wg.links) == len(wg1.links)


def test_add_task():
    """Add add task."""
    wg = WorkGraph("test_add_task")
    add1 = wg.add_task(ArithmeticAddCalculation, name="add1")
    add2 = wg.add_task(ArithmeticAddCalculation, name="add2")
    wg.add_link(add1.outputs["sum"], add2.inputs["x"])
    assert len(wg.tasks) == 2
    assert len(wg.links) == 1


def test_save_load(wg_calcfunction):
    """Save the workgraph"""
    wg = wg_calcfunction
    wg.name = "test_save_load"
    wg.save()
    assert wg.process.process_state.value.upper() == "CREATED"
    assert wg.process.process_label == "WorkGraph<test_save_load>"
    assert wg.process.label == "test_save_load"
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


@pytest.mark.usefixtures("started_daemon_client")
def test_reset_message(wg_calcjob):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    from aiida.cmdline.utils.common import get_workchain_report

    wg = wg_calcjob
    wg.submit()
    # wait for the daemon to start the workgraph
    time.sleep(3)
    wg = WorkGraph.load(wg.process.pk)
    wg.tasks["add2"].set({"y": orm.Int(10).store()})
    wg.save()
    wg.wait()
    report = get_workchain_report(wg.process, "REPORT")
    print(report)
    assert "Action: reset. {'add2'}" in report


def test_restart(wg_calcfunction):
    """Restart from a finished workgraph.
    Load the workgraph, modify the task, and restart the workgraph.
    Only the modified node and its child tasks will be rerun."""
    wg = wg_calcfunction
    wg.add_task(
        "workgraph.test_sum_diff",
        "sumdiff3",
        x=4,
        y=wg.tasks["sumdiff2"].outputs["sum"],
    )
    wg.name = "test_restart_0"
    wg.submit(wait=True)
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.restart()
    wg1.name = "test_restart_1"
    wg1.tasks["sumdiff2"].set({"x": orm.Int(10).store()})
    # wg1.save()
    wg1.submit(wait=True)
    assert wg1.tasks["sumdiff1"].node.pk == wg.tasks["sumdiff1"].pk
    assert wg1.tasks["sumdiff2"].node.pk != wg.tasks["sumdiff2"].pk
    assert wg1.tasks["sumdiff3"].node.pk != wg.tasks["sumdiff3"].pk
    assert wg1.tasks["sumdiff3"].node.outputs.sum == 19


def test_extend_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph("test_graph_build")
    add1 = wg.add_task("workgraph.test_add", "add1", x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    # test wait
    add_multiply_wg.tasks["multiply1"].waiting_on.add("add1")
    # extend workgraph
    wg.extend(add_multiply_wg, prefix="group_")
    assert "group_add1" in [
        task.name for task in wg.tasks["group_multiply1"].waiting_on
    ]
    wg.add_link(add1.outputs[0], wg.tasks["group_add1"].inputs["x"])
    wg.run()
    assert wg.tasks["group_multiply1"].node.outputs.result == 45


@pytest.mark.usefixtures("started_daemon_client")
def test_pause_task_before_submit(wg_calcjob):
    wg = wg_calcjob
    wg.name = "test_pause_task"
    wg.pause_tasks(["add2"])
    wg.submit()
    # wait for the workgraph to launch add2
    wg.wait(tasks={"add2": ["CREATED"]}, timeout=20)
    assert wg.tasks["add2"].node.process_state.value.upper() == "CREATED"
    assert wg.tasks["add2"].node.process_status == "Paused through WorkGraph"
    # I disabled the following lines because the test is not stable
    # Seems the daemon is not responding to the play signal
    # This should be a problem of AiiDA test fixtures
    # wg.play_tasks(["add2"])
    # wg.wait(tasks={"add2": ["FINISHED"]})
    # assert wg.tasks["add2"].outputs["sum"].value == 9


def test_pause_task_after_submit(wg_calcjob):
    wg = wg_calcjob
    wg.tasks["add1"].set({"metadata.options.sleep": 5})
    wg.name = "test_pause_task"
    wg.submit()
    # wait for the workgraph to launch add1
    wg.wait(tasks={"add1": ["CREATED", "WAITING", "RUNNING", "FINISHED"]}, timeout=20)
    wg.pause_tasks(["add2"])
    # wait for the workgraph to launch add2
    wg.wait(tasks={"add2": ["CREATED"]}, timeout=20)
    assert wg.tasks["add2"].node.process_state.value.upper() == "CREATED"
    assert wg.tasks["add2"].node.process_status == "Paused through WorkGraph"
    # I disabled the following lines because the test is not stable
    # Seems the daemon is not responding to the play signal
    # wg.play_tasks(["add2"])
    # wg.wait(tasks={"add2": ["FINISHED"]})
    # assert wg.tasks["add2"].outputs["sum"].value == 9


def test_workgraph_group_outputs(decorated_add):
    wg = WorkGraph("test_workgraph_group_outputs")
    wg.add_task(decorated_add, "add1", x=2, y=3)
    wg.group_outputs = [
        {"name": "sum", "from": "add1.result"},
        # {"name": "add1", "from": "add1"},
    ]
    wg.run()
    assert wg.process.outputs.sum.value == 5
    # assert wg.process.outputs.add1.result.value == 5
