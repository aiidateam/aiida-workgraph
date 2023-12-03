from node_graph.property import NodeProperty as GraphNodeProperty


class NodeProperty(GraphNodeProperty):
    """Represent a property of a Node in the AiiDA WorkTree."""

    @classmethod
    def new(cls, identifier, name=None, data={}):
        """Create a property from a identifier."""
        # use property_pool from aiida_worktree.properties
        # to override the default property_pool from node_graph
        from aiida_worktree.properties import property_pool

        return super().new(
            identifier, name=name, data=data, property_pool=property_pool
        )
