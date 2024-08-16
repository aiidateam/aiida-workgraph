from typing import Dict
from aiida_workgraph.task import Task


class TimeMonitor(Task):
    """Monitor the time"""

    identifier = "workgraph.time_monitor"
    name = "TimeMonitor"
    node_type = "MONITOR"
    catalog = "Control"
    args = ["datetime"]
    kwargs = ["interval"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "datetime")
        inp = self.inputs.new("workgraph.any", "interval")
        inp.add_property("workgraph.any", default=1.0)
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "_wait")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtins",
            "name": "time_monitor",
        }


class FileMonitor(Task):
    """Monitor the file"""

    identifier = "workgraph.file_monitor"
    name = "FileMonitor"
    node_type = "MONITOR"
    catalog = "Control"
    args = ["filepath"]
    kwargs = ["interval"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "filepath")
        inp = self.inputs.new("workgraph.any", "interval")
        inp.add_property("workgraph.any", default=1.0)
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "_wait")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtins",
            "name": "file_monitor",
        }
