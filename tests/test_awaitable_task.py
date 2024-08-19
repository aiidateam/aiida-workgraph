from aiida_workgraph import WorkGraph, task
import asyncio
from aiida.cmdline.utils.common import get_workchain_report
import datetime
import pytest


def test_awaitable_task(decorated_add):
    @task.awaitable()
    async def awaitable_func(x, y):
        n = 2
        while n > 0:
            n -= 1
            await asyncio.sleep(0.5)
        return x + y

    wg = WorkGraph(name="test_awaitable")
    awaitable_func1 = wg.add_task(awaitable_func, "awaitable_func1", x=1, y=2)
    add1 = wg.add_task(decorated_add, "add1", x=1, y=awaitable_func1.outputs["result"])
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: awaitable_func1" in report
    assert add1.outputs["result"].value == 4


def test_time_monitor(decorated_add):
    """Test the time monitor task."""
    wg = WorkGraph(name="test_time_monitor")
    monitor1 = wg.add_task(
        "workgraph.time_monitor",
        "monitor1",
        datetime=datetime.datetime.now() + datetime.timedelta(seconds=10),
    )
    add1 = wg.add_task(decorated_add, "add1", x=1, y=2)
    add1.waiting_on.add(monitor1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: monitor1" in report
    assert add1.outputs["result"].value == 3


def test_file_monitor(decorated_add, tmp_path):
    """Test the file monitor task."""

    @task.awaitable()
    async def create_test_file(filepath="/tmp/test_file_monitor.txt", t=2):
        await asyncio.sleep(t)
        with open(filepath, "w") as f:
            f.write("test")

    monitor_file_path = tmp_path / "test_file_monitor.txt"

    wg = WorkGraph(name="test_file_monitor")
    monitor1 = wg.add_task(
        "workgraph.file_monitor", name="monitor1", filepath=monitor_file_path
    )
    wg.add_task(create_test_file, "create_test_file1", filepath=monitor_file_path)
    add1 = wg.add_task(decorated_add, "add1", x=1, y=2)
    add1.waiting_on.add(monitor1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: monitor1" in report
    assert add1.outputs["result"].value == 3


@pytest.mark.usefixtures("started_daemon_client")
def test_task_monitor(decorated_add):
    """Test the file monitor task."""
    wg2 = WorkGraph(name="wg2")
    monitor1 = wg2.add_task(
        "workgraph.task_monitor",
        name="monitor1",
        workgraph_name="wg1",
        task_name="add1",
    )
    add1 = wg2.add_task(decorated_add, "add1", x=1, y=2, t=0)
    add1.waiting_on.add(monitor1)
    wg2.submit()
    #
    wg1 = WorkGraph(name="wg1")
    wg1.add_task(decorated_add, "add1", x=1, y=2, t=5)
    wg1.submit(wait=True)
    wg2.wait()
    assert wg2.tasks["add1"].node.ctime > wg1.tasks["add1"].node.ctime
