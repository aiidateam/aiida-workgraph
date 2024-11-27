import pytest
from aiida_workgraph import WorkGraph
from typing import Callable
from aiida.cmdline.utils.common import get_workchain_report


def test_task_collection(decorated_add: Callable) -> None:
    """Test the TaskCollection class.
    Since waiting_on and children are TaskCollection instances, we test the waiting_on."""
    wg = WorkGraph("test_task_collection")
    for i in range(5):
        wg.add_task(decorated_add, name=f"task{i}")
    task1 = wg.tasks["task1"]
    # check the graph is not None
    assert task1.waiting_on is not None
    # add a task to waiting_on
    task1.waiting_on.add("task2")
    task1.waiting_on.add(["task3", wg.tasks["task4"]])
    assert len(task1.waiting_on) == 3
    assert wg.tasks["task4"] in task1.waiting_on
    # remove a task from waiting_on
    task1.waiting_on.remove("task2")
    assert wg.tasks["task2"] not in task1.waiting_on
    # clear waiting_on
    task1.waiting_on.clear()
    assert len(task1.waiting_on) == 0


@pytest.mark.usefixtures("started_daemon_client")
def test_build_task_from_workgraph(
    wg_calcfunction: Callable, decorated_add: Callable
) -> None:

    wg = WorkGraph("build_task_from_workgraph")
    add1_task = wg.add_task(decorated_add, name="add1", x=1, y=3)
    wg_task = wg.add_task(wg_calcfunction, name="wg_calcfunction")
    assert wg_task.inputs["sumdiff1"].value is None
    wg.add_task(decorated_add, name="add2", y=3)
    wg.add_link(add1_task.outputs["result"], wg_task.inputs["sumdiff1.x"])
    wg.add_link(wg_task.outputs["sumdiff2.sum"], wg.tasks["add2"].inputs["x"])
    assert len(wg_task.inputs) == 7
    assert len(wg_task.outputs) == 8
    wg.submit(wait=True)
    assert wg.tasks["add2"].outputs["result"].value.value == 14


@pytest.mark.usefixtures("started_daemon_client")
def test_task_wait(decorated_add: Callable) -> None:
    """Run a WorkGraph with a task that waits on other tasks."""

    wg = WorkGraph(name="test_task_wait")
    add1 = wg.add_task(decorated_add, "add1", x=1, y=1)
    add2 = wg.add_task(decorated_add, "add2", x=2, y=2)
    add2.waiting_on.add(add1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add1" in report


def test_set_inputs(decorated_add: Callable) -> None:
    """Test setting inputs of a task."""

    wg = WorkGraph(name="test_set_inputs")
    add1 = wg.add_task(decorated_add, "add1", x=1)
    add1.set({"y": 2, "metadata.store_provenance": False})
    data = wg.prepare_inputs(metadata=None)
    assert data["wg"]["tasks"]["add1"]["inputs"]["y"]["property"]["value"] == 2
    assert (
        data["wg"]["tasks"]["add1"]["inputs"]["metadata"]["property"]["value"][
            "store_provenance"
        ]
        is False
    )
