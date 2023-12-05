from node_graph.node import Node as GraphNode
from aiida_worktree.properties import property_pool
from aiida_worktree.sockets import socket_pool


class Node(GraphNode):
    """Represent a Node in the AiiDA WorkTree.

    The class extends from node_graph.node.Node and add new
    attributes to it.
    """

    property_pool = property_pool
    socket_pool = socket_pool

    def __init__(self, **kwargs):
        """
        Initialize a Node instance.
        """
        super().__init__(**kwargs)
        self.to_ctx = None
        self.wait = None
        self.process = None
        self.pk = None

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

    @classmethod
    def new(cls, identifier, name=None):
        """Create a node from a identifier."""
        from aiida_worktree.nodes import node_pool

        return super().new(identifier, name=name, node_pool=node_pool)

    @classmethod
    def from_dict(cls, data, node_pool=None):
        """Create a node from a dictionary."""
        from aiida_worktree.nodes import node_pool

        node = super().from_dict(data, node_pool=node_pool)
        node.to_ctx = data.get("to_ctx", [])
        node.wait = data.get("wait", [])
        node.process = data.get("process", None)

        return node

    def reset(self):
        self.process = None
        self.state = "CREATED"
