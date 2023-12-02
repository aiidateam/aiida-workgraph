from node_graph.socket import NodeSocket as GraphNodeSocket
from aiida_worktree.property import NodeProperty


class NodeSocket(GraphNodeSocket):
    """Represent a socket of a Node in the AiiDA WorkTree."""

    node_property = NodeProperty
