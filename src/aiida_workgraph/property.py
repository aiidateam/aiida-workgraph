from typing import Union, Optional, Callable
from node_graph.property import NodeProperty


class TaskProperty(NodeProperty):
    """Represent a property of a Task in the AiiDA WorkGraph."""

    def validate(self, value: any) -> None:
        super().validate(value)

    @classmethod
    def new(cls, identifier: Union[Callable, str], name: Optional[str] = None, **kwargs) -> 'TaskProperty':
        """Create a property from a identifier."""
        # use PropertyPool from aiida_workgraph.properties
        # to override the default PropertyPool from node_graph
        from aiida_workgraph.properties import PropertyPool

        return super().new(identifier, name=name, PropertyPool=PropertyPool, **kwargs)
