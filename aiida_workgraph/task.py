from __future__ import annotations

from node_graph.node import Node as GraphNode
from aiida_workgraph import USE_WIDGET
from aiida_workgraph.properties import property_pool
from aiida_workgraph.sockets import socket_pool

if USE_WIDGET:
    from aiida_workgraph.widget import NodeGraphWidget
from aiida_workgraph.collection import (
    WorkGraphPropertyCollection,
    WorkGraphInputSocketCollection,
    WorkGraphOutputSocketCollection,
)
import aiida
from typing import Any, Dict, Optional, Union, Callable, List, Set, Iterable
from aiida_workgraph.utils.message import WIDGET_INSTALLATION_MESSAGE


class Task(GraphNode):
    """Represent a Task in the AiiDA WorkGraph.

    The class extends from node_graph.node.Node and add new
    attributes to it.
    """

    property_pool = property_pool
    socket_pool = socket_pool
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
            input_collection_class=WorkGraphInputSocketCollection,
            output_collection_class=WorkGraphOutputSocketCollection,
            **kwargs,
        )
        self.context_mapping = {} if context_mapping is None else context_mapping
        self.waiting_on = TaskCollection(parent=self)
        self.process = process
        self.pk = pk
        if USE_WIDGET:
            self._widget = NodeGraphWidget(
                settings={"minmap": False},
                style={"width": "80%", "height": "600px"},
            )
        else:
            self._widget = None
        self.state = "PLANNED"
        self.action = ""
        self.show_socket_depth = 0

    def to_dict(self, short: bool = False) -> Dict[str, Any]:
        from aiida.orm.utils.serialize import serialize

        tdata = super().to_dict(short=short)
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
        remain_keys = set(context.values()).difference(self.outputs.keys())
        if remain_keys:
            msg = f"Keys {remain_keys} are not in the outputs of this task."
            raise ValueError(msg)
        self.context_mapping.update(context)

    def set_from_protocol(self, *args: Any, **kwargs: Any) -> None:
        """Set the task inputs from protocol data."""
        from aiida_workgraph.utils import get_executor, get_dict_from_builder

        executor = get_executor(self.get_executor())[0]
        builder = executor.get_builder_from_protocol(*args, **kwargs)
        data = get_dict_from_builder(builder)
        self.set(data)

    @classmethod
    def new(
        cls, identifier: Union[str, Callable], name: Optional[str] = None
    ) -> "Task":
        """Create a task from a identifier."""
        from aiida_workgraph.tasks import task_pool

        return super().new(identifier, name=name, node_pool=task_pool)

    @classmethod
    def from_dict(cls, data: Dict[str, Any], task_pool: Optional[Any] = None) -> "Task":
        """Create a task from a dictionary. This method initializes a Node instance with properties and settings
        defined within the provided data dictionary. If task_pool is not specified, the default task_pool from
        aiida_workgraph.tasks is used.

        Args:
            data (Dict[str, Any]): A dictionary containing the task's configuration.
            task_pool (Optional[Any]): A pool of node configurations, defaults to None
            which will use the global task_pool.

        Returns:
            Node: An instance of Node initialized with the provided data."""
        from aiida_workgraph.tasks import task_pool as workgraph_task_pool
        from aiida.orm.utils.serialize import deserialize_unsafe

        if task_pool is None:
            task_pool = workgraph_task_pool
        task = GraphNode.from_dict(data, node_pool=task_pool)
        task.context_mapping = data.get("context_mapping", {})
        task.waiting_on.add(data.get("wait", []))
        process = data.get("process", None)
        if process and isinstance(process, str):
            process = deserialize_unsafe(process)
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
        from aiida_workgraph.utils import build_callable
        from aiida.engine import ExitCode

        if self._error_handlers is None:
            return {}

        handlers = {}
        if isinstance(self._error_handlers, dict):
            for handler in self._error_handlers.values():
                handler["handler"] = build_callable(handler["handler"])
        elif isinstance(self._error_handlers, list):
            for handler in self._error_handlers:
                handler["handler"] = build_callable(handler["handler"])
                handlers[handler["handler"]["name"]] = handler
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

    def _repr_mimebundle_(self, *args: Any, **kwargs: Any) -> any:

        if self._widget is None:
            print(WIDGET_INSTALLATION_MESSAGE)
            return
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.from_node(self, show_socket_depth=self.show_socket_depth)
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
        if self._widget is None:
            print(WIDGET_INSTALLATION_MESSAGE)
            return
        self._widget.from_node(node=self, show_socket_depth=show_socket_depth)
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
                if task not in self.graph.tasks.keys():
                    raise ValueError(
                        f"Task '{task}' is not in the graph. Available tasks: {self.graph.tasks.keys()}"
                    )
                task_objects.append(self.graph.tasks[task])
            elif isinstance(task, Task):
                task_objects.append(task)
            else:
                raise ValueError(f"Invalid task type: {type(task)}")
        return task_objects

    def add(self, tasks: Union[List[Union[str, Task]], str, Task]) -> None:
        """Add tasks to the collection. Tasks can be a list or a single Task or task name."""
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
