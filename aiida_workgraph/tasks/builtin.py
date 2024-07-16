from typing import Dict
from aiida_workgraph.task import Task


class AiiDAGather(Task):
    """AiiDAGather"""

    identifier = "AiiDAGather"
    name = "AiiDAGather"
    node_type = "WORKCHAIN"
    catalog = "AiiDA"
    kwargs = ["datas"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("Any", "datas")
        inp.link_limit = 100000
        self.outputs.new("Any", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtin",
            "name": "GatherWorkChain",
        }


class AiiDAToCtx(Task):
    """AiiDAToCtx"""

    identifier = "ToCtx"
    name = "ToCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key", "value"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Any", "key")
        self.inputs.new("Any", "value")
        self.outputs.new("Any", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "builtins",
            "name": "setattr",
        }


class AiiDAFromCtx(Task):
    """AiiDAFromCtx"""

    identifier = "FromCtx"
    name = "FromCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Any", "key")
        self.outputs.new("Any", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "builtins",
            "name": "getattr",
        }
