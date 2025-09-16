from .workgraph import WorkGraph
from .task import Task
from .decorator import task
from .tasks import TaskPool
from .tasks.shelljob_task import shelljob
from .utils.flow_control import if_, while_, map_
from .manager import get_current_graph, If, Map, While, Zone
from . import socket_spec as spec
from .socket_spec import namespace, dynamic
from .collection import group

__version__ = "0.7.1"

__all__ = [
    "WorkGraph",
    "Task",
    "task",
    "if_",
    "while_",
    "map_",
    "get_current_graph",
    "Zone",
    "If",
    "Map",
    "While",
    "TaskPool",
    "shelljob",
    "spec",
    "namespace",
    "dynamic",
    "group",
]
