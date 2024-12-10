from aiida_workgraph.socket import TaskSocket


class SocketAny(TaskSocket):
    """Any socket."""

    _socket_identifier: str = "workgraph.any"
    _socket_property_identifier: str = "workgraph.any"


class SocketFloat(TaskSocket):
    """Float socket."""

    _socket_identifier: str = "workgraph.float"
    _socket_property_identifier: str = "workgraph.float"


class SocketInt(TaskSocket):
    """Int socket."""

    _socket_identifier: str = "workgraph.int"
    _socket_property_identifier: str = "workgraph.int"


class SocketString(TaskSocket):
    """String socket."""

    _socket_identifier: str = "workgraph.string"
    _socket_property_identifier: str = "workgraph.string"


class SocketBool(TaskSocket):
    """Bool socket."""

    _socket_identifier: str = "workgraph.bool"
    _socket_property_identifier: str = "workgraph.bool"


class SocketAiiDAFloat(TaskSocket):
    """AiiDAFloat socket."""

    _socket_identifier: str = "workgraph.aiida_float"
    _socket_property_identifier: str = "workgraph.aiida_float"


class SocketAiiDAInt(TaskSocket):
    """AiiDAInt socket."""

    _socket_identifier: str = "workgraph.aiida_int"
    _socket_property_identifier: str = "workgraph.aiida_int"


class SocketAiiDAString(TaskSocket):
    """AiiDAString socket."""

    _socket_identifier: str = "workgraph.aiida_string"
    _socket_property_identifier: str = "workgraph.aiida_string"


class SocketAiiDABool(TaskSocket):
    """AiiDABool socket."""

    _socket_identifier: str = "workgraph.aiida_bool"
    _socket_property_identifier: str = "workgraph.aiida_bool"


class SocketAiiDAList(TaskSocket):
    """AiiDAList socket."""

    _socket_identifier: str = "workgraph.aiida_list"
    _socket_property_identifier: str = "workgraph.aiida_list"


class SocketAiiDADict(TaskSocket):
    """AiiDADict socket."""

    _socket_identifier: str = "workgraph.aiida_dict"
    _socket_property_identifier: str = "workgraph.aiida_dict"


class SocketAiiDAIntVector(TaskSocket):
    """Socket with a AiiDAIntVector property."""

    _socket_identifier: str = "workgraph.aiida_int_vector"
    _socket_property_identifier: str = "workgraph.aiida_int_vector"


class SocketAiiDAFloatVector(TaskSocket):
    """Socket with a FloatVector property."""

    _socket_identifier: str = "workgraph.aiida_float_vector"
    _socket_property_identifier: str = "workgraph.aiida_float_vector"


class SocketStructureData(TaskSocket):
    """Any socket."""

    _socket_identifier: str = "workgraph.aiida_structuredata"
    _socket_property_identifier: str = "workgraph.aiida_structuredata"
