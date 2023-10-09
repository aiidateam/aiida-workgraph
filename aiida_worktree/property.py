from node_graph.property import NodeProperty as GraphNodeProperty


class NodeProperty(GraphNodeProperty):
    """Represent a property of a Node in the AiiDA WorkTree."""

    property_entry = "aiida_worktree.property"
