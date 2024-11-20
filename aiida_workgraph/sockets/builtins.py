from typing import Optional, Any
from aiida_workgraph.socket import TaskSocket


class SocketAny(TaskSocket):
    """Any socket."""

    identifier: str = "workgraph.any"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.any", name, **kwargs)


class SocketNamespace(TaskSocket):
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
        # Set the default value to an empty dictionary
        kwargs.setdefault("default", {})
        self.add_property("workgraph.any", name, **kwargs)


class SocketFloat(TaskSocket):
    """Float socket."""

    identifier: str = "workgraph.float"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.float", name, **kwargs)


class SocketInt(TaskSocket):
    """Int socket."""

    identifier: str = "workgraph.int"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.int", name, **kwargs)


class SocketString(TaskSocket):
    """String socket."""

    identifier: str = "workgraph.string"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.string", name, **kwargs)


class SocketBool(TaskSocket):
    """Bool socket."""

    identifier: str = "workgraph.bool"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.bool", name, **kwargs)


class SocketAiiDAFloat(TaskSocket):
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


class SocketAiiDAInt(TaskSocket):
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


class SocketAiiDAString(TaskSocket):
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


class SocketAiiDABool(TaskSocket):
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


class SocketAiiDAIntVector(TaskSocket):
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


class SocketAiiDAFloatVector(TaskSocket):
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


class SocketStructureData(TaskSocket):
    """Any socket."""

    identifier: str = "workgraph.aiida_structuredata"

    def __init__(
        self, name, node=None, type="INPUT", index=0, uuid=None, **kwargs
    ) -> None:
        super().__init__(name, node, type, index, uuid=uuid)
        self.add_property("workgraph.aiida_structuredata", name, **kwargs)
