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
        from aiida_workgraph.tasks.factory.pythonjob_task import PythonJobTaskFactory
        from aiida_workgraph.tasks.factory.shelljob_task import ShellJobTaskFactory
        from aiida_workgraph.tasks.factory.workgraph_task import WorkGraphTaskFactory

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_task_from_callable(identifier)
        if isinstance(identifier, str) and identifier.upper() == "PYTHONJOB":
            identifier = PythonJobTaskFactory.from_function(kwargs.pop("function"))
        if isinstance(identifier, str) and identifier.upper() == "SHELLJOB":
            TaskCls = ShellJobTaskFactory.create_task(
                outputs=kwargs.get("outputs", None),
                parser_outputs=kwargs.pop("parser_outputs", None),
            )
            task = super()._new(TaskCls, name, uuid, **kwargs)
            return task
        if isinstance(identifier, str) and identifier.upper() == "WHILE":
            task = super()._new("workgraph.while", name, uuid, **kwargs)
            return task
        if isinstance(identifier, str) and identifier.upper() == "IF":
            task = super()._new("workgraph.if", name, uuid, **kwargs)
            return task
        if isinstance(identifier, WorkGraph):
            identifier = WorkGraphTaskFactory.create_task(identifier)
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
