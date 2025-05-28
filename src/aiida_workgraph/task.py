from __future__ import annotations

from node_graph.node import Node as GraphNode
from node_graph.executor import NodeExecutor
from aiida_workgraph.properties import PropertyPool
from aiida_workgraph.sockets import SocketPool
from aiida_workgraph.socket import TaskSocketNamespace
from aiida_workgraph.collection import (
    WorkGraphPropertyCollection,
)
import aiida
from typing import Any, Dict, Optional, Union, Callable, List, Set, Iterable


class Task(GraphNode):
    """Represent a Task in the AiiDA WorkGraph.

    The class extends from node_graph.node.Node and add new
    attributes to it.
    """

    PropertyPool = PropertyPool
    SocketPool = SocketPool
    InputCollectionClass = TaskSocketNamespace
    OutputCollectionClass = TaskSocketNamespace
    PropertyCollectionClass = WorkGraphPropertyCollection

    identifier: str = "workgraph.task"
    is_aiida_component = False
    _error_handlers = None

    def __init__(
        self,
        process: Optional[aiida.orm.ProcessNode] = None,
        pk: Optional[int] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize a Task instance.
        """
        super().__init__(
            **kwargs,
        )
        self.waiting_on = WaitingTaskSet(parent=self)
        self.process = process
        self.pk = pk
        self.state = "PLANNED"
        self.action = ""
        self.show_socket_depth = 0
        self.parent_task = None
        self.map_data = None
        self.mapped_tasks = None
        self.execution_count = 0

    def to_dict(
        self, short: bool = False, should_serialize: bool = False
    ) -> Dict[str, Any]:
        from aiida.orm.utils.serialize import serialize

        tdata = super().to_dict(short=short, should_serialize=should_serialize)
        # clear unused keys
        for key in ["ctrl_inputs", "ctrl_outputs"]:
            tdata.pop(key, None)
        tdata["wait"] = [task.name for task in self.waiting_on]
        tdata["children"] = []
        tdata["execution_count"] = self.execution_count
        tdata["parent_task"] = [self.parent_task.name] if self.parent_task else [None]
        tdata["process"] = serialize(self.process) if self.process else serialize(None)
        tdata["metadata"]["pk"] = self.process.pk if self.process else None
        tdata["metadata"]["is_aiida_component"] = self.is_aiida_component
        tdata["error_handlers"] = self.get_error_handlers()

        return tdata

    def set_from_builder(self, builder: Any) -> None:
        """Set the task inputs from a AiiDA ProcessBuilder."""
        from aiida_workgraph.utils import get_dict_from_builder

        data = get_dict_from_builder(builder)
        self.set(data)

    def set_from_protocol(self, *args: Any, **kwargs: Any) -> None:
        """Set the task inputs from protocol data."""

        executor = NodeExecutor(**self.get_executor()).executor
        # check if the executor has the get_builder_from_protocol method
        if not hasattr(executor, "get_builder_from_protocol"):
            raise AttributeError(
                f"Executor {executor.__name__} does not have the get_builder_from_protocol method."
            )
        builder = executor.get_builder_from_protocol(*args, **kwargs)
        self.set_from_builder(builder)

    @classmethod
    def new(
        cls, identifier: Union[str, Callable], name: Optional[str] = None
    ) -> "Task":
        """Create a task from a identifier."""
        from aiida_workgraph.tasks import TaskPool

        return super().new(identifier, name=name, NodePool=TaskPool)

    @classmethod
    def from_dict(cls, data: Dict[str, Any], TaskPool: Optional[Any] = None) -> "Task":
        """Create a task from a dictionary. This method initializes a Node instance with properties and settings
        defined within the provided data dictionary. If TaskPool is not specified, the default TaskPool from
        aiida_workgraph.tasks is used.

        Args:
            data (Dict[str, Any]): A dictionary containing the task's configuration.
            TaskPool (Optional[Any]): A pool of node configurations, defaults to None
            which will use the global TaskPool.

        Returns:
            Node: An instance of Node initialized with the provided data."""
        from aiida_workgraph.tasks import TaskPool as workgraph_TaskPool

        if TaskPool is None:
            TaskPool = workgraph_TaskPool
        task = GraphNode.from_dict(data, NodePool=TaskPool)

        return task

    def update_from_dict(self, data: Dict[str, Any]) -> None:
        from aiida_workgraph.orm.utils import deserialize_safe

        super().update_from_dict(data)
        process = data.get("process", None)
        if process and isinstance(process, str):
            process = deserialize_safe(process)
        self.process = process
        self._error_handlers = data.get("error_handlers", [])
        self.waiting_on.add(data.get("wait", []))
        self.map_data = data.get("map_data", None)

    def reset(self) -> None:
        self.process = None
        self.state = "PLANNED"

    @property
    def error_handlers(self) -> list:
        return self.get_error_handlers()

    def get_error_handlers(self) -> list:
        """Get the error handler function for this task."""
        from aiida.engine import ExitCode

        if self._error_handlers is None:
            return {}

        handlers = {}
        if isinstance(self._error_handlers, dict):
            for handler in self._error_handlers.values():
                handler["handler"] = NodeExecutor.from_callable(
                    handler["handler"]
                ).to_dict()
        elif isinstance(self._error_handlers, list):
            for handler in self._error_handlers:
                if isinstance(handler["handler"], dict):
                    name = handler.get("name", handler["handler"]["callable_name"])
                else:
                    name = handler.get("name", handler["handler"].__name__)
                    handler["handler"] = NodeExecutor.from_callable(
                        handler["handler"]
                    ).to_dict()
                handlers[name] = handler
        # convert exit code label (str) to status (int)
        for handler in handlers.values():
            exit_codes = []
            for exit_code in handler["exit_codes"]:
                if isinstance(exit_code, int):
                    exit_codes.append(exit_code)
                elif isinstance(exit_code, ExitCode):
                    exit_codes.append(exit_code.status)
                else:
                    raise ValueError(f"Exit code {exit_code} is not a valid exit code.")
            handler["exit_codes"] = exit_codes
        return handlers

    def update_state(self, data: Dict[str, Any]) -> None:
        """Set the outputs of the task from a dictionary."""
        self.state = data["state"]
        self.ctime = data["ctime"]
        self.mtime = data["mtime"]
        self.pk = data["pk"]
        if data["pk"] is not None:
            node = aiida.orm.load_node(data["pk"])
            self.process = self.node = node
            if isinstance(node, aiida.orm.ProcessNode):
                self.set_outputs_from_process_node(node)
            elif isinstance(node, aiida.orm.Data):
                self.set_outputs_from_data_node(node)

    def set_outputs_from_process_node(self, node: aiida.orm.ProcessNode) -> None:
        from aiida_workgraph.utils import get_nested_dict

        # if the node is finished ok, update the output sockets
        # note the task.state may not be the same as the node.process_state
        # for example, task.state can be `SKIPPED` if it is inside a conditional block,
        # even if the node.is_finished_ok is True
        self.process = node
        if node.is_finished_ok:
            # update the output sockets
            for socket in self.outputs:
                if socket._identifier == "workgraph.namespace":
                    socket._value = get_nested_dict(
                        node.outputs, socket._name, default=None
                    )
                else:
                    socket.value = get_nested_dict(
                        node.outputs, socket._name, default=None
                    )

    def set_outputs_from_data_node(self, node: aiida.orm.Data) -> None:
        self.outputs[0].value = node

    def execute(self, args=None, kwargs=None, var_kwargs=None):
        """Execute the task."""
        from node_graph.executor import NodeExecutor

        executor = NodeExecutor(**self.get_executor()).executor
        # the imported executor could be a wrapped function
        if hasattr(executor, "_NodeCls") and hasattr(executor, "_func"):
            executor = getattr(executor, "_func")
        if var_kwargs is None:
            result = executor(*args, **kwargs)
        else:
            result = executor(*args, **kwargs, **var_kwargs)
        return result, "FINISHED"

    def to_widget_value(self):
        from aiida_workgraph.utils import filter_keys_namespace_depth

        tdata = self.to_dict()

        # Remove certain elements of the dict-representation of the Node that we don't want to show
        for key in ("properties", "executor", "node_class", "process"):
            tdata.pop(key, None)

        for input_data in tdata["inputs"]["sockets"].values():
            input_data.pop("property", None)

        tdata["label"] = tdata["identifier"]

        filtered_inputs = filter_keys_namespace_depth(
            dict_=tdata["inputs"]["sockets"], max_depth=self.show_socket_depth
        )
        tdata["inputs"] = list(filtered_inputs.values())
        tdata["outputs"] = list(tdata["outputs"]["sockets"].values())
        wgdata = {"name": self.name, "nodes": {self.name: tdata}, "links": []}
        return wgdata

    def _repr_mimebundle_(self, *args: Any, **kwargs: Any) -> any:
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self.widget.value = self.to_widget_value()
        if hasattr(self.widget, "_repr_mimebundle_"):
            return self.widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self.widget._ipython_display_(*args, **kwargs)

    def to_html(
        self, output: str = None, show_socket_depth: Optional[int] = None, **kwargs
    ):
        """Write a standalone html file to visualize the task."""
        if show_socket_depth is None:
            show_socket_depth = self.show_socket_depth
        self.widget.value = self.to_widget_value()
        return self.widget.to_html(output=output, **kwargs)


class TaskSet:
    """
    A collection of tasks that belong to a specific parent task.

    This class provides an interface for managing a set of tasks, allowing
    tasks to be added, removed, and iterated over while maintaining a reference
    to their parent task and ensuring they exist within a valid task graph.

    """

    def __init__(self, parent: "Task"):
        self._items: Set[str] = set()
        self.parent = parent

    @property
    def graph(self) -> "WorkGraph":
        """Cache and return the graph of the parent task."""
        return self.parent.graph

    @property
    def items(self) -> Set[str]:
        return self._items

    def _normalize_tasks(
        self, tasks: Union[List[Union[str, Task]], str, Task]
    ) -> Iterable[str]:
        """Normalize input to an iterable of task names."""
        if isinstance(tasks, (str, Task)):
            tasks = [tasks]
        task_objects = []
        for task in tasks:
            if isinstance(task, str):
                if task not in self.graph.tasks:
                    raise ValueError(
                        f"Task '{task}' is not in the graph. Available tasks: {self.graph.tasks}"
                    )
                task_objects.append(self.graph.tasks[task])
            elif isinstance(task, Task):
                task_objects.append(task)
            else:
                raise ValueError(f"Invalid task type: {type(task)}")
        return task_objects

    def add(self, tasks: Union[List[Union[str, Task]], str, Task]) -> None:
        """Add tasks to the collection. Tasks can be a list or a single Task or task name."""
        # If the task does not belong to any graph, skip adding it
        if isinstance(self.graph, Task):
            return
        normalize_tasks = []
        for task in self._normalize_tasks(tasks):
            self._items.add(task)
            normalize_tasks.append(task)
        return normalize_tasks

    def remove(self, tasks: Union[List[Union[str, Task]], str, Task]) -> None:
        """Remove tasks from the collection. Tasks can be a list or a single Task or task name."""
        for task in self._normalize_tasks(tasks):
            if task not in self._items:
                raise ValueError(f"Task '{task.name}' is not in the collection.")
            self._items.remove(task)

    def clear(self) -> None:
        """Clear all items from the collection."""
        self._items.clear()

    def __contains__(self, item: str) -> bool:
        """Check if a task name is in the collection."""
        return item in self._items

    def __iter__(self):
        """Yield each task name in the collection for iteration."""
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __repr__(self) -> str:
        return f"{self._items}"


# collection of child tasks
class ChildTaskSet(TaskSet):
    def add(self, tasks: Union[List[Union[str, Task]], str, Task]) -> None:
        """Add tasks to the collection. Tasks can be a list or a single Task or task name."""
        normalize_tasks = super().add(tasks)
        for task in normalize_tasks:
            if task.parent_task is not None:
                raise ValueError(
                    "Task is already a child of the task: {task.parent_task}"
                )
            task.parent_task = self.parent


class WaitingTaskSet(TaskSet):
    def add(self, tasks: Union[List[Union[str, Task]], str, Task]) -> None:
        """Add tasks to the collection. Tasks can be a list or a single Task or task name."""
        normalize_tasks = super().add(tasks)
        for task in normalize_tasks:
            source = task.outputs._wait
            target = self.parent.inputs._wait
            self.graph.add_link(source, target)
