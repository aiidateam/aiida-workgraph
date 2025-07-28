from .workgraph import WorkGraph
from .task import Task
from .decorator import task, build_task
from .tasks import TaskPool
from .tasks.factory.shelljob_task import shelljob
from .utils.flow_control import if_, while_, map_
from .manager import get_current_graph, If, Map, While, Zone

__version__ = "0.6.0"

__all__ = [
    "WorkGraph",
    "Task",
    "task",
    "build_task",
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
]
