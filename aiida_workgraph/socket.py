from typing import Any, Type
from node_graph.socket import NodeSocket
from aiida_workgraph.property import TaskProperty


class TaskSocket(NodeSocket):
    """Represent a socket of a Task in the AiiDA WorkGraph."""

    # use TaskProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    node_property = TaskProperty


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

        def get_serialize(self) -> dict:
            serialize = {"path": "aiida.orm.utils.serialize", "name": "serialize"}
            return serialize

        def get_deserialize(self) -> dict:
            deserialize = {
                "path": "aiida.orm.utils.serialize",
                "name": "deserialize_unsafe",
            }
            return deserialize

    return AiiDATaskSocket
