from typing import Any, Dict
from aiida_workgraph.task import Task, ChildTaskSet


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
        self.children = ChildTaskSet(parent=self)

    def add_task(self, *args, **kwargs):
        """Syntactic sugar to add a task to the zone."""
        task = self.graph.add_task(*args, **kwargs)
        self.children.add(task)
        task.parent_task = self
        return task

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "_wait")

    def to_dict(self, **kwargs) -> Dict[str, Any]:
        tdata = super().to_dict(**kwargs)
        tdata["children"] = [task.name for task in self.children]
        return tdata


class While(Zone):
    """While"""

    identifier = "workgraph.while_zone"
    name = "While"
    node_type = "WHILE"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_input("workgraph.int", "max_iterations", property={"default": 10000})
        self.add_input("workgraph.any", "conditions", link_limit=100000)
        self.add_output("workgraph.any", "_wait")


class If(Zone):
    """If task"""

    identifier = "workgraph.if_zone"
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
            "workgraph.bool", "invert_condition", property={"default": False}
        )
        self.add_output("workgraph.any", "_wait")


class Map(Zone):
    """Map"""

    identifier = "workgraph.map_zone"
    name = "Map"
    node_type = "MAP"
    catalog = "Control"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.graph is not None:
            # init a ctx variable to placeholder the item
            key = f"map_zone_{self.name}_item"
            self.graph.update_ctx({key: None})
            self.item = self.graph.ctx[key]
        else:
            raise ValueError("Map zone must be added to a graph.")

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_input("workgraph.any", "source", link_limit=100000)
        self.add_input("workgraph.any", "placeholder")
        self.add_output("workgraph.any", "_wait")


class SetContext(Task):
    """SetContext"""

    identifier = "workgraph.set_context"
    name = "SetContext"
    node_type = "Normal"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "context")
        self.add_input("workgraph.any", "key")
        self.add_input("workgraph.any", "value")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.builtins",
            "callable_name": "update_ctx",
        }
        return executor


class GetContext(Task):
    """GetContext"""

    identifier = "workgraph.get_context"
    name = "GetContext"
    node_type = "Normal"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "context")
        self.add_input("workgraph.any", "key")
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.builtins",
            "callable_name": "get_context",
        }
        return executor


class AiiDAInt(Task):
    identifier = "workgraph.aiida_int"
    name = "AiiDAInt"
    node_type = "Normal"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.add_input("workgraph.any", "value", property={"default": 0.0})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.int", "result")
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
    node_type = "Normal"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.float", "value", property={"default": 0.0})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.float", "result")
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
    node_type = "Normal"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.string", "value", property={"default": ""})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.string", "result")
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
    node_type = "Normal"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "value", property={"default": []})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.list", "result")
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
    node_type = "Normal"
    catalog = "Test"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "value", property={"default": {}})
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.dict", "result")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida.orm",
            "callable_name": "Dict",
        }
        return executor


class AiiDANode(Task):
    """AiiDANode"""

    identifier = "workgraph.load_node"
    name = "AiiDANode"
    node_type = "Normal"
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

    identifier = "workgraph.load_code"
    name = "AiiDACode"
    node_type = "Normal"
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


class GraphBuilderTask(Task):
    """Graph builder task"""

    identifier = "workgraph.graph_builder"
    name = "graph_builder"
    node_type = "graph_builder"
    catalog = "builtins"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine
        from node_graph.executor import NodeExecutor

        executor = NodeExecutor(**self.get_executor()).executor

        if var_kwargs is None:
            wg = executor(*args, **kwargs)
        else:
            wg = executor(*args, **kwargs, **var_kwargs)
        wg.name = self.name

        wg.parent_uuid = engine_process.node.uuid
        inputs = wg.prepare_inputs(metadata={"call_link_label": self.name})
        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                WorkGraphEngine,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(WorkGraphEngine, **inputs)
            state = "RUNNING"
        process.label = self.name

        return process, state
