from aiida_workgraph.socket import TaskSocket


class SocketAny(TaskSocket):
    """Any socket."""

    identifier: str = "workgraph.any"
    property_identifier: str = "workgraph.any"


class SocketNamespace(TaskSocket):
    """Namespace socket."""

    identifier: str = "workgraph.namespace"
    property_identifier: str = "workgraph.any"


class SocketFloat(TaskSocket):
    """Float socket."""

    identifier: str = "workgraph.float"
    property_identifier: str = "workgraph.float"


class SocketInt(TaskSocket):
    """Int socket."""

    identifier: str = "workgraph.int"
    property_identifier: str = "workgraph.int"


class SocketString(TaskSocket):
    """String socket."""

    identifier: str = "workgraph.string"
    property_identifier: str = "workgraph.string"


class SocketBool(TaskSocket):
    """Bool socket."""

    identifier: str = "workgraph.bool"
    property_identifier: str = "workgraph.bool"


class SocketAiiDAFloat(TaskSocket):
    """AiiDAFloat socket."""

    identifier: str = "workgraph.aiida_float"
    property_identifier: str = "workgraph.aiida_float"


class SocketAiiDAInt(TaskSocket):
    """AiiDAInt socket."""

    identifier: str = "workgraph.aiida_int"
    property_identifier: str = "workgraph.aiida_int"


class SocketAiiDAString(TaskSocket):
    """AiiDAString socket."""

    identifier: str = "workgraph.aiida_string"
    property_identifier: str = "workgraph.aiida_string"


class SocketAiiDABool(TaskSocket):
    """AiiDABool socket."""

    identifier: str = "workgraph.aiida_bool"
    property_identifier: str = "workgraph.aiida_bool"


class SocketAiiDAIntVector(TaskSocket):
    """Socket with a AiiDAIntVector property."""

    identifier: str = "workgraph.aiida_int_vector"
    property_identifier: str = "workgraph.aiida_int_vector"


class SocketAiiDAFloatVector(TaskSocket):
    """Socket with a FloatVector property."""

    identifier: str = "workgraph.aiida_float_vector"
    property_identifier: str = "workgraph.aiida_float_vector"


class SocketStructureData(TaskSocket):
    """Any socket."""

    identifier: str = "workgraph.aiida_structuredata"
    property_identifier: str = "workgraph.aiida_structuredata"
