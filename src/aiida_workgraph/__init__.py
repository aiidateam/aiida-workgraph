from .workgraph import WorkGraph
from .task import Task
from .decorator import task, build_task
from .tasks import TaskPool
from .utils.flow_control import if_, while_, map_
from .manager import active_graph, If, Map, While

__version__ = "0.5.3"

__all__ = [
    "WorkGraph",
    "Task",
    "task",
    "build_task",
    "if_",
    "while_",
    "map_",
    "active_graph",
    "If",
    "Map",
    "While",
    "TaskPool",
]
