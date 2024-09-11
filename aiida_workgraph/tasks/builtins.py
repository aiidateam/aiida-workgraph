from typing import Any, Dict
from aiida_workgraph.task import Task, TaskCollection


class Zone(Task):
    """
    Extend the Task class to include a 'children' attribute.
    """

    identifier = "workgraph.zone"
    name = "Zone"
    node_type = "ZONE"
    catalog = "Control"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.children = TaskCollection(parent=self)

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "_wait")

    def to_dict(self) -> Dict[str, Any]:
        tdata = super().to_dict()
        tdata["children"] = [task.name for task in self.children]
        return tdata

    def from_dict(self, data: Dict[str, Any]) -> None:
        super().from_dict(data)
        self.children.add(data.get("children", []))


class While(Zone):
    """While"""

    identifier = "workgraph.while"
    name = "While"
    node_type = "WHILE"
    catalog = "Control"
    kwargs = ["max_iterations", "conditions"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        inp = self.inputs.new("node_graph.int", "max_iterations")
        inp.add_property("node_graph.int", default=10000)
        inp = self.inputs.new("workgraph.any", "conditions")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "_wait")


class If(Zone):
    """If task"""

    identifier = "workgraph.if"
    name = "If"
    node_type = "IF"
    catalog = "Control"
    kwargs = ["conditions", "invert_condition"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.inputs.new("workgraph.any", "conditions")
        self.inputs.new("workgraph.any", "invert_condition")
        self.outputs.new("workgraph.any", "_wait")


class Gather(Task):
    """Gather"""

    identifier = "workgraph.aiida_gather"
    name = "Gather"
    node_type = "WORKCHAIN"
    catalog = "Control"
    kwargs = ["datas"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("workgraph.any", "datas")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtins",
            "name": "GatherWorkChain",
        }


class ToContext(Task):
    """ToContext"""

    identifier = "workgraph.to_context"
    name = "ToContext"
    node_type = "TO_CONTEXT"
    catalog = "Control"
    args = ["key", "value"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        self.inputs.new("workgraph.any", "value")


class FromContext(Task):
    """FromContext"""

    identifier = "workgraph.from_context"
    name = "FromContext"
    node_type = "FROM_CONTEXT"
    catalog = "Control"
    args = ["key"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        self.outputs.new("workgraph.any", "result")


class AiiDAInt(Task):
    identifier = "workgraph.aiida_int"
    name = "AiiDAInt"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_sockets(self) -> None:
        inp = self.inputs.new("workgraph.any", "value", default=0.0)
        inp.add_property("workgraph.aiida_int", default=1.0)
        self.outputs.new("workgraph.aiida_int", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Int",
        }


class AiiDAFloat(Task):
    identifier = "workgraph.aiida_float"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_sockets(self) -> None:
        self.inputs.new("workgraph.aiida_float", "value", default=0.0)
        self.outputs.new("workgraph.aiida_float", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Float",
        }


class AiiDAString(Task):
    identifier = "workgraph.aiida_string"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_sockets(self) -> None:
        self.inputs.new("AiiDAString", "value", default="")
        self.outputs.new("AiiDAString", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Str",
        }


class AiiDAList(Task):
    identifier = "workgraph.aiida_list"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseList", "value", default=[])

    def create_sockets(self) -> None:
        self.outputs.new("workgraph.any", "Parameters")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "List",
        }


class AiiDADict(Task):
    identifier = "workgraph.aiida_dict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseDict", "value", default={})

    def create_sockets(self) -> None:
        self.outputs.new("workgraph.any", "Parameters")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Dict",
        }


class AiiDANode(Task):
    """AiiDANode"""

    identifier = "workgraph.aiida_node"
    name = "AiiDANode"
    node_type = "node"
    catalog = "Test"
    kwargs = ["identifier", "pk", "uuid", "label"]

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "identifier")
        self.inputs.new("workgraph.any", "pk")
        self.inputs.new("workgraph.any", "uuid")
        self.inputs.new("workgraph.any", "label")
        self.outputs.new("workgraph.any", "node")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "load_node",
        }


class AiiDACode(Task):
    """AiiDACode"""

    identifier = "workgraph.aiida_code"
    name = "AiiDACode"
    node_type = "node"
    catalog = "Test"
    kwargs = ["identifier", "pk", "uuid", "label"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "identifier")
        self.inputs.new("workgraph.any", "pk")
        self.inputs.new("workgraph.any", "uuid")
        self.inputs.new("workgraph.any", "label")
        self.outputs.new("workgraph.any", "Code")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "load_code",
        }


class Select(Task):
    """Select"""

    identifier = "workgraph.select"
    name = "Select"
    node_type = "Normal"
    catalog = "Control"
    args = ["condition", "true", "false"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "condition")
        self.inputs.new("workgraph.any", "true")
        self.inputs.new("workgraph.any", "false")
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "result")
        self.outputs.new("workgraph.any", "_wait")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtins",
            "name": "select",
        }
