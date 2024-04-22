from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
    InputSocketCollection,
    OutputSocketCollection,
)


class WorkGraphNodeCollection(NodeCollection):
    def new(self, identifier, name=None, uuid=None, **kwargs):
        from aiida_workgraph.decorator import build_node_from_callable

        # build the node on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_node_from_callable(identifier)
        # Call the original new method
        return super().new(identifier, name, uuid, **kwargs)


class WorkGraphPropertyCollection(PropertyCollection):
    def new(self, identifier, name=None, **kwargs):
        from aiida_workgraph.property import build_property_from_AiiDA

        # build the property on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_property_from_AiiDA(identifier)
        # Call the original new method
        return super().new(identifier, name, **kwargs)


class WorkGraphInputSocketCollection(InputSocketCollection):
    def new(self, identifier, name=None, **kwargs):
        from aiida_workgraph.socket import build_socket_from_AiiDA

        # build the socket on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_socket_from_AiiDA(identifier)
        # Call the original new method
        return super().new(identifier, name, **kwargs)


class WorkGraphOutputSocketCollection(OutputSocketCollection):
    def new(self, identifier, name=None, **kwargs):
        from aiida_workgraph.socket import build_socket_from_AiiDA

        # build the socket on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_socket_from_AiiDA(identifier)
        # Call the original new method
        return super().new(identifier, name, **kwargs)
