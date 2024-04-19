from node_graph.socket import NodeSocket as GraphNodeSocket
from aiida_workgraph.property import NodeProperty


class NodeSocket(GraphNodeSocket):
    """Represent a socket of a Node in the AiiDA WorkGraph."""

    # use NodeProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    node_property = NodeProperty
