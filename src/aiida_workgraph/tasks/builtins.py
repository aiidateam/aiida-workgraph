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
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
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

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.inputs.new(
            "node_graph.int", "max_iterations", property_data={"default": 10000}
        )
        self.inputs.new("workgraph.any", "conditions", link_limit=100000)
        self.outputs.new("workgraph.any", "_wait")


class If(Zone):
    """If task"""

    identifier = "workgraph.if"
    name = "If"
    node_type = "IF"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.inputs.new("workgraph.any", "conditions")
        self.inputs.new("workgraph.any", "invert_condition")
        self.outputs.new("workgraph.any", "_wait")


class SetContext(Task):
    """SetContext"""

    identifier = "workgraph.set_context"
    name = "SetContext"
    node_type = "SET_CONTEXT"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        self.inputs.new("workgraph.any", "value")
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.any", "_wait")


class GetContext(Task):
    """GetContext"""

    identifier = "workgraph.get_context"
    name = "GetContext"
    node_type = "GET_CONTEXT"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "key")
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
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

    def create_sockets(self) -> None:
        self.inputs.new("workgraph.any", "value", property_data={"default": 0.0})
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.aiida_int", "result")
        self.outputs.new("workgraph.any", "_wait")


class AiiDAFloat(Task):
    identifier = "workgraph.aiida_float"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Float",
    }

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.float", "value", property_data={"default": 0.0})
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.aiida_float", "result")
        self.outputs.new("workgraph.any", "_wait")


class AiiDAString(Task):
    identifier = "workgraph.aiida_string"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Str",
    }

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.string", "value", property_data={"default": ""})
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.aiida_string", "result")
        self.outputs.new("workgraph.any", "_wait")


class AiiDAList(Task):
    identifier = "workgraph.aiida_list"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "List",
    }

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "value", property_data={"default": []})
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.aiida_list", "result")
        self.outputs.new("workgraph.any", "_wait")


class AiiDADict(Task):
    identifier = "workgraph.aiida_dict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    _executor = {
        "module": "aiida.orm",
        "name": "Dict",
    }

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "value", property_data={"default": {}})
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.aiida_dict", "result")
        self.outputs.new("workgraph.any", "_wait")


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

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "identifier")
        self.inputs.new("workgraph.any", "pk")
        self.inputs.new("workgraph.any", "uuid")
        self.inputs.new("workgraph.any", "label")
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.any", "node")
        self.outputs.new("workgraph.any", "_wait")


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

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "identifier")
        self.inputs.new("workgraph.any", "pk")
        self.inputs.new("workgraph.any", "uuid")
        self.inputs.new("workgraph.any", "label")
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.any", "Code")
        self.outputs.new("workgraph.any", "_wait")


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

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("workgraph.any", "condition")
        self.inputs.new("workgraph.any", "true")
        self.inputs.new("workgraph.any", "false")
        self.inputs.new("workgraph.any", "_wait", arg_type="none", link_limit=100000)
        self.outputs.new("workgraph.any", "result")
        self.outputs.new("workgraph.any", "_wait")
