try:
    import anywidget  # noqa: F401

    USE_WIDGET = True
except ImportError:
    USE_WIDGET = False

from .workgraph import WorkGraph
from .task import Task
from .decorator import task, build_task

__version__ = "0.3.14"

__all__ = ["WorkGraph", "Task", "task", "build_task"]
