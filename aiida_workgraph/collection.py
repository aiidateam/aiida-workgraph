from node_graph.collection import NodeCollection


class WorkGraphNodeCollection(NodeCollection):
    def new(self, identifier, name=None, uuid=None, **kwargs):
        from aiida_workgraph.decorator import build_node_from_callable

        # build the node on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_node_from_callable(identifier)
        # Call the original new method
        return super().new(identifier, name, uuid, **kwargs)
