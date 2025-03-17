"""
Simple global variable approach.
Note pitfalls:
    - lack of concurrency control
    - collisions in library code
    - difficulty testing in isolation
    - etc.
"""
from contextlib import contextmanager


class CurrentGraphManager:
    _instance = None

    def __new__(cls, *args, **kwargs):
        """
        Enforce the singleton pattern. Only one instance of
        CurrentGraphManager is created for the entire process.
        """
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._graph = None  # Storage for the active graph
        return cls._instance

    def get_current_graph(self):
        """
        Retrieve the current graph, or create a new one if none is set.
        """
        from aiida_workgraph import WorkGraph

        if self._graph is None:
            self._graph = WorkGraph()
        return self._graph

    def set_current_graph(self, graph):
        """
        Set the active graph to the given instance.
        """
        self._graph = graph

    @contextmanager
    def active_graph(self, graph):
        """
        Context manager that temporarily overrides the current graph
        with `graph`, restoring the old graph when exiting the context.
        """
        old_graph = self._graph
        self._graph = graph
        try:
            yield graph
        finally:
            self._graph = old_graph


# Create a global manager instance
_manager = CurrentGraphManager()


def get_current_graph():
    """
    Helper function to retrieve the graph
    through the global manager instance.
    """
    return _manager.get_current_graph()


def set_current_graph(graph):
    """
    Helper function to set the graph through the
    global manager instance.
    """
    _manager.set_current_graph(graph)


@contextmanager
def active_graph(graph):
    """
    Top-level context manager that defers to
    the manager's `active_graph` method.
    """
    with _manager.active_graph(graph) as g:
        yield g
