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
        self.to_ctx = []
        self.wait = []
        self.process = None

    def to_dict(self):
        ndata = super().to_dict()
        ndata["to_ctx"] = self.to_ctx
        ndata["wait"] = self.wait
        ndata["process"] = self.process.uuid if self.process else None

        return ndata
