from typing import Any, Dict
from aiida_workgraph.task import Task, ChildTaskSet
from aiida_workgraph.tasks.factory.base import BaseTaskFactory


class GraphLevelTask(Task):
    """Base class for graph level tasks"""

    catalog = "Builtins"
    is_dynamic: bool = True
    node_class = Task
    factory_class = BaseTaskFactory

    @property
    def outputs(self):
        return self.inputs

    @outputs.setter
    def outputs(self, _value):
        """Outputs are the same as inputs for ctx node."""
        pass

    def get_metadata(self):

        metadata = super().get_metadata()
        metadata["node_class"] = {
            "module_path": self.node_class.__module__,
            "callable_name": self.node_class.__name__,
        }
        metadata["factory_class"] = {
            "module_path": self.factory_class.__module__,
            "callable_name": self.factory_class.__name__,
        }
        return metadata


class GraphInputs(GraphLevelTask):
    identifier = "workgraph.graph_inputs"
    name = "Graph_Inputs"


class GraphOutputs(GraphLevelTask):
    identifier = "workgraph.graph_outputs"
    name = "Graph_Outputs"


class GraphContext(GraphLevelTask):
    identifier = "workgraph.graph_ctx"
    name = "Graph_Ctx"


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
        task.parent = self
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

    @property
    def item(self):
        for child in self.children:
            if child.identifier == "workgraph.map_item":
                return child.outputs.item
        # create a child map_item_task if it does not exist
        map_item_task = self.add_task("workgraph.map_item")
        return map_item_task.outputs.item

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "source", link_limit=100000)
        self.add_input(
            "workgraph.any", "_wait", link_limit=100000, metadata={"arg_type": "none"}
        )
        self.add_output("workgraph.any", "_wait")


class MapItem(Task):
    """MapItem"""

    identifier = "workgraph.map_item"
    name = "MapItem"
    node_type = "Normal"
    catalog = "Control"

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "source", link_limit=100000)
        self.add_input("workgraph.any", "key")
        self.add_output("workgraph.any", "item")
        self.add_output("workgraph.any", "_wait")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.builtins",
            "callable_name": "get_item",
        }
        return executor


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


class GraphTask(Task):
    """Graph builder task"""

    identifier = "workgraph.graph_task"
    name = "graph_task"
    node_type = "graph_task"
    catalog = "builtins"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine
        from aiida_workgraph.decorator import _run_func_with_wg
        from node_graph.executor import NodeExecutor
        from aiida_workgraph import task

        executor = NodeExecutor(**self.get_executor()).executor
        # Cloudpickle doesn’t restore the function’s own name in its globals after unpickling,
        # so any recursive calls would raise NameError. As a temporary workaround, we re-insert
        # the decorated function into its globals under its original name.
        # Downside: this mutates the module globals at runtime, if another symbol with the same name exists,
        # we may introduce hard-to-trace bugs or collisions.
        if executor.__name__ not in executor.__globals__:
            if getattr(executor, "is_decoratored", False):
                executor.__globals__[executor.__name__] = executor
            else:
                executor.__globals__[executor.__name__] = task.graph()(executor)
        if getattr(executor, "is_decoratored", False):
            executor = executor._func
        graph_task_output_names = [
            output._name
            for output in self.outputs
            if not output._metadata.builtin_socket
        ]
        wg = _run_func_with_wg(
            executor, graph_task_output_names, args, kwargs, var_kwargs
        )
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
