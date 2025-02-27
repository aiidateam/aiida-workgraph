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

    def add_task(self, *args, **kwargs):
        """Syntactic sugar to add a task to the zone."""
        task = self.parent.add_task(*args, **kwargs)
        self.children.add(task)
        return task

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "_wait")

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
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_input(
            "node_graph.int", "max_iterations", property_data={"default": 10000}
        )
        self.add_input("workgraph.any", "conditions", link_limit=100000)
        self.add_output("workgraph.any", "_wait")


class If(Zone):
    """If task"""

    identifier = "workgraph.if"
    name = "If"
    node_type = "IF"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_input("workgraph.any", "conditions")
        self.add_input(
            "workgraph.bool", "invert_condition", property_data={"default": False}
        )
        self.add_output("workgraph.any", "_wait")


class SetContext(Task):
    """SetContext"""

    identifier = "workgraph.set_context"
    name = "SetContext"
    node_type = "SET_CONTEXT"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "key")
        self.add_input("workgraph.any", "value")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "_wait")


class GetContext(Task):
    """GetContext"""

    identifier = "workgraph.get_context"
    name = "GetContext"
    node_type = "GET_CONTEXT"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "key")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")


class AiiDAInt(Task):
    identifier = "workgraph.aiida_int"
    name = "AiiDAInt"
    node_type = "data"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.add_input("workgraph.any", "value", property_data={"default": 0.0})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.aiida_int", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "Int",
        }
        return executor


class AiiDAFloat(Task):
    identifier = "workgraph.aiida_float"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.float", "value", property_data={"default": 0.0})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.aiida_float", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "Float",
        }
        return executor


class AiiDAString(Task):
    identifier = "workgraph.aiida_string"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.string", "value", property_data={"default": ""})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.aiida_string", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "Str",
        }
        return executor


class AiiDAList(Task):
    identifier = "workgraph.aiida_list"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "value", property_data={"default": []})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.aiida_list", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "List",
        }
        return executor


class AiiDADict(Task):
    identifier = "workgraph.aiida_dict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "value", property_data={"default": {}})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.aiida_dict", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "Dict",
        }
        return executor


class AiiDANode(Task):
    """AiiDANode"""

    identifier = "workgraph.aiida_node"
    name = "AiiDANode"
    node_type = "node"
    catalog = "Test"

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "identifier")
        self.add_input("workgraph.any", "pk")
        self.add_input("workgraph.any", "uuid")
        self.add_input("workgraph.any", "label")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "node")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "load_node",
        }
        return executor


class AiiDACode(Task):
    """AiiDACode"""

    identifier = "workgraph.aiida_code"
    name = "AiiDACode"
    node_type = "node"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "identifier")
        self.add_input("workgraph.any", "pk")
        self.add_input("workgraph.any", "uuid")
        self.add_input("workgraph.any", "label")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "Code")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "load_code",
        }
        return executor


class Select(Task):
    """Select"""

    identifier = "workgraph.select"
    name = "Select"
    node_type = "Normal"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "condition")
        self.add_input("workgraph.any", "true")
        self.add_input("workgraph.any", "false")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.builtins",
            "callable_name": "select",
        }
        return executor


class WorkGraphTask(Task):

    identifier = "workgraph.workgraph"
    name = "AiiDAWorkGraph"
    node_type = "workgraph"
    catalog = "WORKGRAPH"
    # is_aiida_component = True  # ???

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._workgraph = None

    @property
    def workgraph(self):
        from aiida_workgraph import WorkGraph
        if not self._workgraph:
            graph_data = self.get_executor()["graph_data"]
            self._workgraph = WorkGraph.from_dict(graph_data)
        return self._workgraph

    @property
    def tasks(self):
        return self.workgraph.tasks

    @property
    def links(self):
        return self.workgraph.links
