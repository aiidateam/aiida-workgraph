from typing import Optional, Any
from aiida_workgraph.socket import TaskSocket
from node_graph.serializer import SerializeJson, SerializePickle


class SocketAny(TaskSocket, SerializePickle):
    """Any socket."""

    identifier: str = "workgraph.any"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.any", name, **kwargs)


class SocketNamespace(TaskSocket, SerializePickle):
    """Namespace socket."""

    identifier: str = "workgraph.namespace"

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
        self.add_property("workgraph.any", name, **kwargs)


class SocketAiiDAFloat(TaskSocket, SerializeJson):
    """AiiDAFloat socket."""

    identifier: str = "workgraph.aiida_float"

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
        self.add_property("workgraph.aiida_float", name, **kwargs)


class SocketAiiDAInt(TaskSocket, SerializeJson):
    """AiiDAInt socket."""

    identifier: str = "workgraph.aiida_int"

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
        self.add_property("workgraph.aiida_int", name, **kwargs)


class SocketAiiDAString(TaskSocket, SerializeJson):
    """AiiDAString socket."""

    identifier: str = "workgraph.aiida_string"

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
        self.add_property("workgraph.aiida_string", name, **kwargs)


class SocketAiiDABool(TaskSocket, SerializeJson):
    """AiiDABool socket."""

    identifier: str = "workgraph.aiida_bool"

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
        self.add_property("workgraph.aiida_bool", name, **kwargs)


class SocketAiiDAIntVector(TaskSocket, SerializeJson):
    """Socket with a AiiDAIntVector property."""

    identifier: str = "workgraph.aiida_int_vector"

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
        self.add_property("workgraph.aiida_int_vector", name, **kwargs)


class SocketAiiDAFloatVector(TaskSocket, SerializeJson):
    """Socket with a FloatVector property."""

    identifier: str = "workgraph.aiida_float_vector"

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
        self.add_property("workgraph.aiida_float_vector", name, **kwargs)
