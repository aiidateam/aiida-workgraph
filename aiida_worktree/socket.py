from node_graph.socket import NodeSocket as GraphNodeSocket


class NodeSocket(GraphNodeSocket):
    """Represent a socket of a Node in the AiiDA WorkTree."""

    socket_entry = "aiida_worktree.socket"
    property_entry = "aiida_worktree.property"
