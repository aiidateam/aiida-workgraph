from node_graph.node import Node as GraphNode
from aiida_workgraph.properties import property_pool
from aiida_workgraph.sockets import socket_pool
from aiida_workgraph.widget import NodeGraphWidget
from aiida_workgraph.collection import (
    WorkGraphPropertyCollection,
    WorkGraphInputSocketCollection,
    WorkGraphOutputSocketCollection,
)
import aiida
from typing import Any, Dict, Optional, Union, Callable, List


class Task(GraphNode):
    """Represent a Task in the AiiDA WorkGraph.

    The class extends from node_graph.node.Node and add new
    attributes to it.
    """

    property_pool = property_pool
    socket_pool = socket_pool
    is_aiida_component = False

    def __init__(
        self,
        to_context: Optional[List[Any]] = None,
        wait: List[Union[str, GraphNode]] = [],
        process: Optional[aiida.orm.ProcessNode] = None,
        pk: Optional[int] = None,
        **kwargs: Any
    ) -> None:
        """
        Initialize a Task instance.
        """
        super().__init__(
            property_collection_class=WorkGraphPropertyCollection,
            input_collection_class=WorkGraphInputSocketCollection,
            output_collection_class=WorkGraphOutputSocketCollection,
            **kwargs
        )
        self.to_context = [] if to_context is None else to_context
        self.wait = [] if wait is None else wait
        self.process = process
        self.pk = pk
        self._widget = NodeGraphWidget(
            settings={"minmap": False},
            style={"width": "80%", "height": "600px"},
        )

    def to_dict(self) -> Dict[str, Any]:
        tdata = super().to_dict()
        tdata["to_context"] = [] if self.to_context is None else self.to_context
        tdata["wait"] = [
            task if isinstance(task, str) else task.name for task in self.wait
        ]
        tdata["process"] = self.process.uuid if self.process else None
        tdata["metadata"]["pk"] = self.process.pk if self.process else None
        tdata["metadata"]["is_aiida_component"] = self.is_aiida_component

        return tdata

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
        from aiida_workgraph.tasks import task_pool

        task = super().from_dict(data, node_pool=task_pool)
        task.to_context = data.get("to_context", [])
        task.wait = data.get("wait", [])
        task.process = data.get("process", None)

        return task

    def reset(self) -> None:
        self.process = None
        self.state = "CREATED"

    def _repr_mimebundle_(self, *args: Any, **kwargs: Any) -> any:
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.from_node(self)
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)

    def to_html(self, output: str = None, **kwargs):
        """Write a standalone html file to visualize the task."""
        self._widget.from_node(self)
        return self._widget.to_html(output=output, **kwargs)
