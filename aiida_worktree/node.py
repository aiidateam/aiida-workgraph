from node_graph.node import Node as GraphNode


class Node(GraphNode):
    """Represent a Node in the AiiDA WorkTree.

    The class extends from node_graph.node.Node and add new
    attributes to it.
    """

    socket_entry = "aiida_worktree.socket"
    property_entry = "aiida_worktree.property"

    def __init__(self, **kwargs):
        """
        Initialize a Node instance.
        """
        super().__init__(**kwargs)
        self.to_ctx = None
        self.wait = None
        self.process = None

    def to_dict(self):
        ndata = super().to_dict()
        ndata["to_ctx"] = [] if self.to_ctx is None else self.to_ctx
        ndata["wait"] = [] if self.wait is None else self.wait
        ndata["process"] = self.process.uuid if self.process else None

        return ndata

    def set_from_protocol(self, *args, **kwargs):
        """For node support protocol, set the node from protocol data."""
        from aiida_worktree.utils import get_executor, get_dict_from_builder

        executor = get_executor(self.get_executor())[0]
        builder = executor.get_builder_from_protocol(*args, **kwargs)
        data = get_dict_from_builder(builder)
        self.set(data)
