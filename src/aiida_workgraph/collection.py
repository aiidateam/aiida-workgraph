from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
)
from typing import Any, Callable, Optional, Union
from aiida.engine import CalcJob, WorkChain
from aiida_workgraph.utils import build_callable


class TaskCollection(NodeCollection):
    def _new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        uuid: Optional[str] = None,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.decorator import (
            build_task_from_callable,
            build_pythonjob_task,
            build_shelljob_task,
            build_task_from_workgraph,
        )
        from aiida_workgraph.workgraph import WorkGraph

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            if isinstance(identifier, type):
                if issubclass(identifier, (CalcJob, WorkChain)):
                    executor = build_callable(identifier)
                    task = super()._new(
                        "workgraph.aiida_process",
                        name,
                        uuid,
                        _executor=executor,
                        **kwargs
                    )
                    return task
            identifier = build_task_from_callable(identifier)
        if isinstance(identifier, str) and identifier.upper() == "PYTHONJOB":
            identifier, _ = build_pythonjob_task(kwargs.pop("function"))
        if isinstance(identifier, str) and identifier.upper() == "SHELLJOB":
            identifier, _ = build_shelljob_task(
                outputs=kwargs.get("outputs", None),
                parser_outputs=kwargs.pop("parser_outputs", None),
            )
            task = super()._new(identifier, name, uuid, **kwargs)
            return task
        if isinstance(identifier, str) and identifier.upper() == "WHILE":
            task = super()._new("workgraph.while", name, uuid, **kwargs)
            return task
        if isinstance(identifier, str) and identifier.upper() == "IF":
            task = super()._new("workgraph.if", name, uuid, **kwargs)
            return task
        if isinstance(identifier, WorkGraph):
            identifier = build_task_from_workgraph(identifier)
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
