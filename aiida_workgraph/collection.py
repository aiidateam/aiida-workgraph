from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
    InputSocketCollection,
    OutputSocketCollection,
)
from typing import Any, Callable, Optional, Union


class TaskCollection(NodeCollection):
    def new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        uuid: Optional[str] = None,
        run_remotely: Optional[bool] = False,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.decorator import (
            build_task_from_callable,
            build_python_task,
            build_shell_task,
        )

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_task_from_callable(identifier)
            if run_remotely:
                if identifier.node.node_type.upper() == "GRAPH_BUILDER":
                    raise ValueError(
                        "GraphBuilder task cannot be run remotely. Please set run_remotely=False."
                    )
                # this is a PythonTask
                identifier, _ = build_python_task(identifier)
            return super().new(identifier, name, uuid, **kwargs)
        if isinstance(identifier, str) and identifier.upper() == "PythonTask":
            identifier, _ = build_python_task(kwargs.pop("function"))
            return super().new(identifier, name, uuid, **kwargs)
        if isinstance(identifier, str) and identifier.upper() == "SHELLTASK":
            identifier, _, links = build_shell_task(
                nodes=kwargs.get("nodes", {}),
                outputs=kwargs.get("outputs", None),
                parser_outputs=kwargs.pop("parser_outputs", None),
            )
            task = super().new(identifier, name, uuid, **kwargs)
            # make links between the tasks
            task.set(links)
            return task
        return super().new(identifier, name, uuid, **kwargs)


class WorkGraphPropertyCollection(PropertyCollection):
    def new(
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
        return super().new(identifier, name, **kwargs)


class WorkGraphInputSocketCollection(InputSocketCollection):
    def new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.socket import build_socket_from_AiiDA

        # build the socket on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_socket_from_AiiDA(identifier)
        # Call the original new method
        return super().new(identifier, name, **kwargs)


class WorkGraphOutputSocketCollection(OutputSocketCollection):
    def new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.socket import build_socket_from_AiiDA

        # build the socket on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_socket_from_AiiDA(identifier)
        # Call the original new method
        return super().new(identifier, name, **kwargs)
