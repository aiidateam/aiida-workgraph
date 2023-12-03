from node_graph.socket import NodeSocket as GraphNodeSocket
from aiida_worktree.property import NodeProperty


class NodeSocket(GraphNodeSocket):
    """Represent a socket of a Node in the AiiDA WorkTree."""

    # use NodeProperty from aiida_worktree.property
    # to override the default NodeProperty from node_graph
    node_property = NodeProperty
