from node_graph.node import Node as GraphNode
from aiida_workgraph.properties import property_pool
from aiida_workgraph.sockets import socket_pool
from aiida_workgraph.widget import NodeGraphWidget


class Node(GraphNode):
    """Represent a Node in the AiiDA WorkGraph.

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
        self.wait = []
        self.process = None
        self.pk = None
        self._widget = NodeGraphWidget(
            settings={"minmap": False}, style={"width": "40%", "height": "600px"}
        )

    def to_dict(self):
        ndata = super().to_dict()
        ndata["to_ctx"] = [] if self.to_ctx is None else self.to_ctx
        ndata["wait"] = [
            node if isinstance(node, str) else node.name for node in self.wait
        ]
        ndata["process"] = self.process.uuid if self.process else None
        ndata["metadata"]["pk"] = self.process.pk if self.process else None

        return ndata

    def set_from_protocol(self, *args, **kwargs):
        """For node support protocol, set the node from protocol data."""
        from aiida_workgraph.utils import get_executor, get_dict_from_builder

        executor = get_executor(self.get_executor())[0]
        builder = executor.get_builder_from_protocol(*args, **kwargs)
        data = get_dict_from_builder(builder)
        self.set(data)

    @classmethod
    def new(cls, identifier, name=None):
        """Create a node from a identifier."""
        from aiida_workgraph.nodes import node_pool

        return super().new(identifier, name=name, node_pool=node_pool)

    @classmethod
    def from_dict(cls, data, node_pool=None):
        """Create a node from a dictionary."""
        from aiida_workgraph.nodes import node_pool

        node = super().from_dict(data, node_pool=node_pool)
        node.to_ctx = data.get("to_ctx", [])
        node.wait = data.get("wait", [])
        node.process = data.get("process", None)

        return node

    def reset(self):
        self.process = None
        self.state = "CREATED"

    def _repr_mimebundle_(self, *args, **kwargs):
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.from_node(self)
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)
