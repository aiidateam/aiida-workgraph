import asyncio
import datetime
from aiida_pythonjob import MonitorPyFunction
import pytest
from aiida_workgraph import WorkGraph, task
from aiida_workgraph.tasks.monitors import monitor_time


def test_monitor_decorator():
    wg = WorkGraph(name='test_monitor_decorator')
    t = datetime.datetime.now() + datetime.timedelta(seconds=5)
    monitor1 = wg.add_task(
        monitor_time,
        'time_monitor1',
        time=t,
    )
    wg.run()
    assert wg.process.is_finished_ok is True
    assert monitor1.node.is_finished_ok is True
    assert (monitor1.node.mtime - monitor1.node.ctime) >= datetime.timedelta(seconds=5)


def test_monitor_timeout_and_interval():
    with WorkGraph() as wg:
        monitor_time(
            time=datetime.datetime.now() + datetime.timedelta(seconds=20),
            interval=2,
            timeout=5,
        )
        wg.run()
        node = wg.tasks[-1].node
        assert not node.is_finished_ok
        assert node.exit_status == MonitorPyFunction.exit_codes.ERROR_TIMEOUT.status


def test_builtin_time_monitor_entrypoint():
    """Test the time monitor task."""
    wg = WorkGraph(name='test_time_monitor')
    monitor1 = wg.add_task(
        'workgraph.monitor_time',
        'monitor1',
        time=datetime.datetime.now() + datetime.timedelta(seconds=5),
    )
    assert monitor1.get_executor().callable == monitor_time


def test_builtin_file_monitor_entrypoint(tmp_path):
    """Test the file monitor task."""

    @task
    async def create_test_file(filepath='/tmp/test_file_monitor.txt', t=2):
        await asyncio.sleep(t)
        with open(filepath, 'w') as f:
            f.write('test')

    monitor_file_path = str(tmp_path / 'test_file_monitor.txt')

    t = 2
    wg = WorkGraph(name='test_file_monitor')
    wg.add_task(create_test_file, 'create_test_file1', filepath=monitor_file_path, t=t)
    monitor1 = wg.add_task('workgraph.monitor_file', name='monitor1', filepath=monitor_file_path)
    wg.run()
    assert wg.process.is_finished_ok is True
    assert (monitor1.node.mtime - monitor1.node.ctime) >= datetime.timedelta(seconds=t)


@pytest.mark.usefixtures('started_daemon_client')
def test_builtin_task_monitor_entrypoint(decorated_add):
    """Test the file monitor task."""
    wg2 = WorkGraph(name='wg2')
    monitor1 = wg2.add_task(
        'workgraph.monitor_task',
        name='monitor1',
        workgraph_name='wg1',
        task_name='add1',
    )
    add1 = wg2.add_task(decorated_add, 'add1', x=1, y=2, t=0)
    monitor1 >> add1
    wg2.submit()
    #
    wg1 = WorkGraph(name='wg1')
    wg1.add_task(decorated_add, 'add1', x=1, y=2, t=5)
    wg1.run()
    wg2.wait()
    assert wg2.tasks.add1.node.ctime > wg1.tasks.add1.node.ctime


@pytest.mark.usefixtures('started_daemon_client')
def test_builtin_task_monitor_entrypoint_timeout(decorated_add):
    """Test the monitor task with a timeout."""
    wg = WorkGraph(name='test_monitor_timeout')
    monitor1 = wg.add_task(
        'workgraph.monitor_file',
        name='monitor1',
        timeout=2,
        filepath='/tmp/test_file_monitor.txt',
    )
    wg.run()
    node = monitor1.node
    assert not node.is_finished_ok
    assert node.exit_status == MonitorPyFunction.exit_codes.ERROR_TIMEOUT.status


@pytest.mark.usefixtures('started_daemon_client')
def test_task_monitor_kill():
    """Test killing a monitor task."""
    wg = WorkGraph(name='test_monitor_kill')
    monitor1 = wg.add_task(
        'workgraph.monitor_file',
        name='monitor1',
        timeout=30,
        filepath='/tmp/test_file_monitor.txt',
    )
    wg.submit()
    wg.wait(tasks={'monitor1': ['RUNNING']})
    wg.kill_tasks(['monitor1'])
    wg.wait()
    assert monitor1.node.is_killed
