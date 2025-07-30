import datetime
import typing as t

from aiida_workgraph import task
from aiida_workgraph.tasks.factory.awaitable_task import MonitorFunctionTask


async def monitor(function, interval, timeout, *args, **kwargs):
    """Monitor the function until it returns `True` or the timeout is reached."""
    import asyncio
    import time

    start_time = time.time()
    while True:
        result = function(*args, **kwargs)
        if result:
            break
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Timeout reached for monitor function {function}")
        await asyncio.sleep(interval)


@task.monitor
def monitor_file(filepath: str):
    """Return `True` when the file is detected."""
    import os

    return os.path.exists(filepath)


@task.monitor
def monitor_time(time: t.Union[str, datetime.datetime]):
    """Return `True` when the given moment in time has passed.

    If given as a string, `time` should be in ISO format (e.g., '2025-07-30T12:00:00').
    """

    if isinstance(time, str):
        try:
            time = datetime.datetime.fromisoformat(time)
        except ValueError as err:
            raise ValueError(
                f"Invalid time format: {time}. Expected ISO format."
            ) from err

    return datetime.datetime.now() > time


@task.monitor
def monitor_task(task_name: str, workgraph_pk: int = None, workgraph_name: str = None):
    """Return `True` if the task in the WorkGraph is completed."""
    from aiida import orm

    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    if workgraph_pk:
        try:
            node = orm.load_node(workgraph_pk)
        except Exception:
            return False
    else:
        builder = orm.QueryBuilder()
        builder.append(
            WorkGraphEngine,
            filters={
                "attributes.process_label": {"==": f"WorkGraph<{workgraph_name}>"}
            },
            tag="process",
        )
        if builder.count() == 0:
            return False
    print("Found workgraph")
    node = builder.first()[0]
    state = node.task_states.get(task_name, "")
    print(f"Task state: {state}")
    if state in ["FINISHED", "FAILED", "SKIPPED"]:
        return True
    else:
        return False


class TimeMonitor(MonitorFunctionTask):
    """Monitor the time"""

    identifier = "workgraph.time_monitor"
    name = "TimeMonitor"
    catalog = "Monitor"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "time")
        inp = self.add_input("workgraph.any", "interval")
        inp.add_property("workgraph.any", default=1.0)
        inp = self.add_input("workgraph.any", "timeout")
        inp.add_property("workgraph.any", default=86400.0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        inp._link_limit = 100000
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.tasks.monitors",
            "callable_name": "monitor_time",
        }
        return executor


class FileMonitor(MonitorFunctionTask):
    """Monitor the file"""

    identifier = "workgraph.file_monitor"
    name = "FileMonitor"
    catalog = "Monitor"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "filepath")
        inp = self.add_input("workgraph.any", "interval")
        inp.add_property("workgraph.any", default=1.0)
        inp = self.add_input("workgraph.any", "timeout")
        inp.add_property("workgraph.any", default=86400.0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        inp._link_limit = 100000
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.tasks.monitors",
            "callable_name": "monitor_file",
        }
        return executor


class TaskMonitor(MonitorFunctionTask):
    """Monitor the file"""

    identifier = "workgraph.task_monitor"
    name = "TaskMonitor"
    catalog = "Monitor"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "workgraph_pk")
        self.add_input("workgraph.any", "workgraph_name")
        self.add_input("workgraph.any", "task_name")
        inp = self.add_input("workgraph.any", "interval")
        inp.add_property("workgraph.any", default=1.0)
        inp = self.add_input("workgraph.any", "timeout")
        inp.add_property("workgraph.any", default=86400.0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        inp._link_limit = 100000
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.tasks.monitors",
            "callable_name": "monitor_task",
        }
        return executor
