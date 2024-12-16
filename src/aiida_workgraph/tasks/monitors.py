from aiida_workgraph.task import Task


class TimeMonitor(Task):
    """Monitor the time"""

    identifier = "workgraph.time_monitor"
    name = "TimeMonitor"
    node_type = "MONITOR"
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
            "module_path": "aiida_workgraph.executors.monitors",
            "callable_name": "time_monitor",
        }
        return executor


class FileMonitor(Task):
    """Monitor the file"""

    identifier = "workgraph.file_monitor"
    name = "FileMonitor"
    node_type = "MONITOR"
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
            "module_path": "aiida_workgraph.executors.monitors",
            "callable_name": "file_monitor",
        }
        return executor


class TaskMonitor(Task):
    """Monitor the file"""

    identifier = "workgraph.task_monitor"
    name = "TaskMonitor"
    node_type = "MONITOR"
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
            "module_path": "aiida_workgraph.executors.monitors",
            "callable_name": "task_monitor",
        }
        return executor
