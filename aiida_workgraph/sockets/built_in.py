from typing import Optional, Any
from aiida_workgraph.socket import TaskSocket
from node_graph.serializer import SerializeJson, SerializePickle
from node_graph.sockets.builtin import (
    SocketInt,
    SocketFloat,
    SocketString,
    SocketBool,
    SocketBaseDict,
    SocketBaseList,
)


class SocketAny(TaskSocket, SerializePickle):
    """Socket for any time."""

    identifier: str = "Any"

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
        self.add_property("Any", name, **kwargs)


class SocketNamespace(TaskSocket, SerializePickle):
    """Namespace socket."""

    identifier: str = "Namespace"

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
        self.add_property("Any", name, **kwargs)


class SocketAiiDAFloat(TaskSocket, SerializeJson):
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


class SocketAiiDAInt(TaskSocket, SerializeJson):
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


class SocketAiiDAString(TaskSocket, SerializeJson):
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


class SocketAiiDABool(TaskSocket, SerializeJson):
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


class SocketAiiDAIntVector(TaskSocket, SerializeJson):
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


class SocketAiiDAFloatVector(TaskSocket, SerializeJson):
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
    SocketAny,
    SocketNamespace,
    SocketInt,
    SocketFloat,
    SocketString,
    SocketBool,
    SocketBaseDict,
    SocketBaseList,
    SocketAiiDAInt,
    SocketAiiDAFloat,
    SocketAiiDAString,
    SocketAiiDABool,
    SocketAiiDAIntVector,
    SocketAiiDAFloatVector,
]
