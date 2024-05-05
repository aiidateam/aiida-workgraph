from typing import Any, Type
from node_graph.socket import NodeSocket as GraphNodeSocket
from aiida_workgraph.property import NodeProperty


class NodeSocket(GraphNodeSocket):
    """Represent a socket of a Node in the AiiDA WorkGraph."""

    # use NodeProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    node_property = NodeProperty


def build_socket_from_AiiDA(DataClass: Type[Any]) -> Type[GraphNodeSocket]:
    """Create a socket class from AiiDA DataClass."""

    class AiiDANodeSocket(NodeSocket):
        """AiiDA Node Socket."""

        identifier: str = DataClass.__name__

        def __init__(
            self,
            name: str,
            node: Any = None,
            type: str = "INPUT",
            index: int = 0,
            uuid: str = None,
            **kwargs: Any
        ) -> None:
            super().__init__(name, node, type, index, uuid=uuid)
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

    return AiiDANodeSocket
