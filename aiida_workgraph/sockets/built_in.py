from typing import Optional, Any
from aiida_workgraph.socket import NodeSocket
from node_graph.serializer import SerializeJson, SerializePickle
from node_graph.sockets.builtin import (
    SocketBaseDict,
    SocketBaseList,
)


class SocketGeneral(NodeSocket, SerializePickle):
    """General socket."""

    identifier: str = "General"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("General", name, **kwargs)


class SocketAiiDAFloat(NodeSocket, SerializeJson):
    """AiiDAFloat socket."""

    identifier: str = "AiiDAFloat"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("AiiDAFloat", name, **kwargs)


class SocketAiiDAInt(NodeSocket, SerializeJson):
    """AiiDAInt socket."""

    identifier: str = "AiiDAInt"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("AiiDAInt", name, **kwargs)


class SocketAiiDAString(NodeSocket, SerializeJson):
    """AiiDAString socket."""

    identifier: str = "AiiDAString"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("AiiDAString", name, **kwargs)


class SocketAiiDABool(NodeSocket, SerializeJson):
    """AiiDABool socket."""

    identifier: str = "AiiDABool"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("AiiDABool", name, **kwargs)


class SocketAiiDAIntVector(NodeSocket, SerializeJson):
    """Socket with a AiiDAIntVector property."""

    identifier: str = "AiiDAIntVector"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("AiiDAIntVector", name, **kwargs)


class SocketAiiDAFloatVector(NodeSocket, SerializeJson):
    """Socket with a FloatVector property."""

    identifier: str = "FloatVector"

    def __init__(
        self,
        name: str,
        node: Optional[Any] = None,
        type: str = "INPUT",
        index: int = 0,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("FloatVector", name, **kwargs)


socket_list = [
    SocketGeneral,
    SocketBaseDict,
    SocketBaseList,
    SocketAiiDAInt,
    SocketAiiDAFloat,
    SocketAiiDAString,
    SocketAiiDABool,
    SocketAiiDAIntVector,
    SocketAiiDAFloatVector,
]
