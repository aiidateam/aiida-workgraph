from typing import Any, Type

from aiida import orm
from node_graph.socket import NodeSocket

from aiida_workgraph.property import TaskProperty


class TaskSocket(NodeSocket):
    """Represent a socket of a Task in the AiiDA WorkGraph."""

    # use TaskProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    node_property = TaskProperty

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
