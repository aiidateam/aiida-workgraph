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

    def to_dict(self, short: bool = False) -> Dict[str, Any]:
        tdata = super().to_dict(short=short)
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


class SetContext(Task):
    """SetContext"""

    identifier = "workgraph.set_context"
    name = "SetContext"
    node_type = "SET_CONTEXT"
    catalog = "Control"
    args = ["key", "value"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.inputs.new("workgraph.any", "value")
        self.outputs.new("workgraph.any", "_wait")


class GetContext(Task):
    """GetContext"""

    identifier = "workgraph.get_context"
    name = "GetContext"
    node_type = "GET_CONTEXT"
    catalog = "Control"
    args = ["key"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        inp = self.inputs.new("workgraph.any", "_wait")
        inp.link_limit = 100000
        self.outputs.new("workgraph.any", "result")
        self.outputs.new("workgraph.any", "_wait")


class AiiDAInt(Task):
    identifier = "workgraph.aiida_int"
    name = "AiiDAInt"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Int",
    }
    args = ["value"]

    def create_sockets(self) -> None:
        inp = self.inputs.new("workgraph.any", "value", default=0.0)
        inp.add_property("workgraph.aiida_int", default=1.0)
        self.outputs.new("workgraph.aiida_int", "result")


class AiiDAFloat(Task):
    identifier = "workgraph.aiida_float"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Float",
    }
    args = ["value"]

    def create_sockets(self) -> None:
        self.inputs.new(
            "workgraph.aiida_float", "value", property_data={"default": 0.0}
        )
        self.outputs.new("workgraph.aiida_float", "result")


class AiiDAString(Task):
    identifier = "workgraph.aiida_string"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Str",
    }
    args = ["value"]

    def create_sockets(self) -> None:
        self.inputs.new("AiiDAString", "value", default="")
        self.outputs.new("AiiDAString", "result")


class AiiDAList(Task):
    identifier = "workgraph.aiida_list"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "List",
    }
    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseList", "value", default=[])

    def create_sockets(self) -> None:
        self.outputs.new("workgraph.any", "Parameters")


class AiiDADict(Task):
    identifier = "workgraph.aiida_dict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Dict",
    }
    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseDict", "value", default={})

    def create_sockets(self) -> None:
        self.outputs.new("workgraph.any", "Parameters")


class AiiDANode(Task):
    """AiiDANode"""

    identifier = "workgraph.aiida_node"
    name = "AiiDANode"
    node_type = "node"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "load_node",
    }
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


class AiiDACode(Task):
    """AiiDACode"""

    identifier = "workgraph.aiida_code"
    name = "AiiDACode"
    node_type = "node"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "load_code",
    }
    kwargs = ["identifier", "pk", "uuid", "label"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "identifier")
        self.inputs.new("workgraph.any", "pk")
        self.inputs.new("workgraph.any", "uuid")
        self.inputs.new("workgraph.any", "label")
        self.outputs.new("workgraph.any", "Code")


class Select(Task):
    """Select"""

    identifier = "workgraph.select"
    name = "Select"
    node_type = "Normal"
    catalog = "Control"

    _executor = {
        "module": "aiida_workgraph.executors.builtins",
        "name": "select",
    }
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
