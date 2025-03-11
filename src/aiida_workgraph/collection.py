from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
)
from typing import Any, Callable, Optional, Union


class TaskCollection(NodeCollection):
    """Collection of tasks."""


class WorkGraphPropertyCollection(PropertyCollection):
    def _new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.property import build_property_from_AiiDA

        # build the property on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_property_from_AiiDA(identifier)
        # Call the original new method
        return super()._new(identifier, name, **kwargs)
