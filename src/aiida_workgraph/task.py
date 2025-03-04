from __future__ import annotations

from node_graph.node import Node as GraphNode
from node_graph.executor import NodeExecutor
from aiida_workgraph.properties import PropertyPool
from aiida_workgraph.sockets import SocketPool
from aiida_workgraph.socket import NodeSocketNamespace
from node_graph_widget import NodeGraphWidget
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
    is_aiida_component = False
    _error_handlers = None

    def __init__(
        self,
        context_mapping: Optional[List[Any]] = None,
        process: Optional[aiida.orm.ProcessNode] = None,
        pk: Optional[int] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize a Task instance.
        """
        super().__init__(
            property_collection_class=WorkGraphPropertyCollection,
            input_collection_class=NodeSocketNamespace,
            output_collection_class=NodeSocketNamespace,
            **kwargs,
        )
        self.context_mapping = {} if context_mapping is None else context_mapping
        self.waiting_on = TaskCollection(parent=self)
        self.process = process
        self.pk = pk
        self._widget = NodeGraphWidget(
            settings={"minmap": False},
            style={"width": "80%", "height": "600px"},
        )
        self.state = "PLANNED"
        self.action = ""
        self.show_socket_depth = 0

    def to_dict(self, short: bool = False) -> Dict[str, Any]:
        from aiida.orm.utils.serialize import serialize

        tdata = super().to_dict(short=short)
        # clear unused keys
        for key in ["ctrl_inputs", "ctrl_outputs"]:
            tdata.pop(key, None)
        tdata["context_mapping"] = self.context_mapping
        tdata["wait"] = [task.name for task in self.waiting_on]
        tdata["children"] = []
        tdata["execution_count"] = 0
        tdata["parent_task"] = [None]
        tdata["process"] = serialize(self.process) if self.process else serialize(None)
        tdata["metadata"]["pk"] = self.process.pk if self.process else None
        tdata["metadata"]["is_aiida_component"] = self.is_aiida_component
        tdata["error_handlers"] = self.get_error_handlers()

        return tdata

    def set_context(self, context: Dict[str, Any]) -> None:
        """Set the output of the task as a value in the context.
        key is the context key, value is the output key.
        """
        # all values should belong to the outputs.keys()
        remain_keys = set(context.values()).difference(self.get_output_names())
        if remain_keys:
            msg = f"Keys {remain_keys} are not in the outputs of this task."
            raise ValueError(msg)
        self.context_mapping.update(context)

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
        from aiida_workgraph.orm.utils import deserialize_safe

        if TaskPool is None:
            TaskPool = workgraph_TaskPool
        task = GraphNode.from_dict(data, NodePool=TaskPool)
        task.context_mapping = data.get("context_mapping", {})
        task.waiting_on.add(data.get("wait", []))
        process = data.get("process", None)
        if process and isinstance(process, str):
            process = deserialize_safe(process)
        task.process = process
        task._error_handlers = data.get("error_handlers", [])

        return task

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

    def to_widget_value(self):
        from aiida_workgraph.utils import filter_keys_namespace_depth

        tdata = self.to_dict()

        # Remove certain elements of the dict-representation of the Node that we don't want to show
        for key in ("properties", "executor", "node_class", "process"):
            tdata.pop(key, None)
        for input in tdata["inputs"].values():
            input.pop("property", None)

        tdata["label"] = tdata["identifier"]

        filtered_inputs = filter_keys_namespace_depth(
            dict_=tdata["inputs"], max_depth=self.show_socket_depth
        )
        tdata["inputs"] = list(filtered_inputs.values())
        tdata["outputs"] = list(tdata["outputs"].values())
        wgdata = {"name": self.name, "nodes": {self.name: tdata}, "links": []}
        return wgdata

    def _repr_mimebundle_(self, *args: Any, **kwargs: Any) -> any:
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.value = self.to_widget_value()
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)

    def to_html(
        self, output: str = None, show_socket_depth: Optional[int] = None, **kwargs
    ):
        """Write a standalone html file to visualize the task."""
        if show_socket_depth is None:
            show_socket_depth = self.show_socket_depth
        self._widget.value = self.to_widget_value()
        return self._widget.to_html(output=output, **kwargs)


class TaskCollection:
    def __init__(self, parent: "Task"):
        self._items: Set[str] = set()
        self.parent = parent
        self._top_parent = None

    @property
    def graph(self) -> "WorkGraph":
        """Cache and return the top parent of the collection."""
        if not self._top_parent:
            parent = self.parent
            while getattr(parent, "parent", None):
                parent = parent.parent
            self._top_parent = parent
        return self._top_parent

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
        for task in self._normalize_tasks(tasks):
            self._items.add(task)

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
