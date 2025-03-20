from typing import Any, Type

from aiida import orm
from node_graph.socket import (
    NodeSocket,
    NodeSocketNamespace,
)

from aiida_workgraph.property import TaskProperty
from aiida_workgraph.orm.mapping import type_mapping


class TaskSocket(NodeSocket):
    """Represent a socket of a Task in the AiiDA WorkGraph."""

    # use TaskProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    _socket_property_class = TaskProperty

    @property
    def _decorator(self):
        from aiida_workgraph.decorator import task

        return task

    @property
    def node_value(self):
        return self.get_node_value()

    def get_node_value(self):
        """Obtain the actual Python `value` of the object attached to the Socket."""
        if isinstance(self.value, orm.Data):
            if hasattr(self.value, "value"):
                return self.value.value
            else:
                raise ValueError(
                    "Data node does not have a value attribute. We do not know how to extract the raw Python value."
                )
        else:
            return self.value


class TaskSocketNamespace(NodeSocketNamespace):
    """Represent a namespace of a Task in the AiiDA WorkGraph."""

    _identifier = "workgraph.namespace"
    _socket_property_class = TaskProperty
    _type_mapping: dict = type_mapping

    @property
    def _decorator(self):
        from aiida_workgraph.decorator import task

        return task

    def __init__(self, *args, **kwargs):
        super().__init__(*args, entry_point="aiida_workgraph.socket", **kwargs)


class ContextSocketNamespace(TaskSocketNamespace):
    """Namespace used for the context"""

    _identifier = "workgraph.context_namespace"

    def __init__(self, *args, **kwargs):
        from aiida_workgraph.task import Task

        super().__init__(*args, **kwargs)
        # Add a temporary node, but not in the WorkGraph
        # we should disable users adding task with name `_context`
        self._node = Task(name="ctx", graph=self._graph)

    def __setattr__(self, name: str, value: Any) -> None:
        """
        Override __setattr__ so that doing `namespace_socket.some_name = x`
        either sets the property or links to another socket, rather than
        replacing the entire sub-socket object.
        """
        # If the attribute name is "private" or reserved, do normal attribute setting
        if name.startswith("_"):
            object.__setattr__(self, name, value)
            return

        self._set_socket_value({name: value}, link_limit=1e6)


def build_socket_from_AiiDA(DataClass: Type[Any]) -> Type[TaskSocket]:
    """Create a socket class from AiiDA DataClass."""

    class AiiDATaskSocket(TaskSocket):
        """AiiDA Task Socket."""

        identifier: str = DataClass.__name__

        def __init__(
            self,
            name: str,
            parent: Any = None,
            type: str = "INPUT",
            index: int = 0,
            uuid: str = None,
            **kwargs: Any
        ) -> None:
            super().__init__(name, parent, type, index, uuid=uuid)
            self.add_property(DataClass, name, **kwargs)

    return AiiDATaskSocket
