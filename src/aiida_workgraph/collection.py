from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
)
from typing import Any, Callable, Optional, Union


class TaskCollection(NodeCollection):
    def _new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.decorator import build_task_from_callable
        from aiida_workgraph.workgraph import WorkGraph
        from aiida_workgraph.tasks.factory.workgraph_task import WorkGraphTaskFactory

        # build the task on the fly if the identifier is a callable
        if isinstance(identifier, WorkGraph):
            identifier = WorkGraphTaskFactory.create_task(identifier)
        elif callable(identifier):
            identifier = build_task_from_callable(identifier)

        return super()._new(identifier, name, uuid, **kwargs)


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
