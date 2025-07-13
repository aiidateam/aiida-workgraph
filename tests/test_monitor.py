from aiida_workgraph import WorkGraph, task
import asyncio
from aiida.cmdline.utils.common import get_workchain_report
import datetime
import pytest
import time


@task.monitor()
def time_monitor(time: str):
    """Return True if the current time is greater than the given time."""
    import datetime

    # load the time string
    time = datetime.datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f")
    print(f"Monitoring time: {time}")
    return datetime.datetime.now() > time


def test_monitor_decorator():

    wg = WorkGraph(name="test_monitor_decorator")
    wg.add_task(
        time_monitor,
        "time_monitor1",
        time=str(datetime.datetime.now() + datetime.timedelta(seconds=5)),
    )
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: time_monitor1" in report
    assert wg.process.is_finished_ok is True


def test_monitor_timeout_and_interval():

    with WorkGraph() as wg:
        time_monitor(
            time=str(datetime.datetime.now() + datetime.timedelta(seconds=20)),
            interval=2,
            timeout=5,
        )
        wg.run()
        report = get_workchain_report(wg.process, "REPORT")
        assert "Timeout reached for monitor" in report
        assert wg.process.is_finished_ok is False


def test_builtin_time_monitor(decorated_add):
    """Test the time monitor task."""
    wg = WorkGraph(name="test_time_monitor")
    monitor1 = wg.add_task(
        "workgraph.time_monitor",
        "monitor1",
        time=str(datetime.datetime.now() + datetime.timedelta(seconds=5)),
    )
    add1 = wg.add_task(decorated_add, "add1", x=1, y=2)
    add1.waiting_on.add(monitor1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: monitor1" in report
    assert add1.outputs.result.value == 3


def test_builtin_file_monitor(decorated_add, tmp_path):
    """Test the file monitor task."""

    @task.awaitable()
    async def create_test_file(filepath="/tmp/test_file_monitor.txt", t=2):
        await asyncio.sleep(t)
        with open(filepath, "w") as f:
            f.write("test")

    monitor_file_path = str(tmp_path / "test_file_monitor.txt")

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
    assert add1.outputs.result.value == 3


@pytest.mark.usefixtures("started_daemon_client")
def test_builtin_task_monitor(decorated_add):
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
    wg1.submit(wait=True, timeout=60)
    wg2.wait()
    assert wg2.tasks.add1.node.ctime > wg1.tasks.add1.node.ctime


@pytest.mark.usefixtures("started_daemon_client")
def test_builtin_task_monitor_timeout(decorated_add):
    """Test the monitor task with a timeout."""
    wg = WorkGraph(name="test_monitor_timeout")
    monitor1 = wg.add_task(
        "workgraph.file_monitor",
        name="monitor1",
        timeout=2,
        filepath="/tmp/test_file_monitor.txt",
    )
    add1 = wg.add_task(decorated_add, "add1", x=1, y=2)
    add1.waiting_on.add(monitor1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Timeout reached for monitor function" in report
    assert monitor1.state == "FAILED"
    assert add1.state == "SKIPPED"


@pytest.mark.usefixtures("started_daemon_client")
def test_task_monitor_kill(decorated_add):
    """Test killing a monitor task."""
    wg = WorkGraph(name="test_monitor_kill")
    monitor1 = wg.add_task(
        "workgraph.file_monitor",
        name="monitor1",
        timeout=30,
        filepath="/tmp/test_file_monitor.txt",
    )
    add1 = wg.add_task(decorated_add, "add1", x=1, y=2)
    add1.waiting_on.add(monitor1)
    wg.submit()
    time.sleep(5)
    wg.wait(tasks={"monitor1": ["RUNNING"]})
    wg.kill_tasks(["monitor1"])
    wg.wait()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Task monitor1 was KILLED" in report
    assert monitor1.state == "KILLED"
    assert add1.state == "SKIPPED"
