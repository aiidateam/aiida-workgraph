from node_graph.collection import (
    NodeCollection,
    PropertyCollection,
    InputSocketCollection,
    OutputSocketCollection,
)
from typing import Any, Callable, Optional, Union


class WorkGraphNodeCollection(NodeCollection):
    def new(
        self,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        uuid: Optional[str] = None,
        on_remote: Optional[bool] = False,
        **kwargs: Any
    ) -> Any:
        from aiida_workgraph.decorator import (
            build_task_from_callable,
            build_PythonJob_task,
            build_ShellJob_task,
        )

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_task_from_callable(identifier)
            if kwargs.pop("run_remotely", False):
                if identifier.node.node_type.upper() == "GRAPH_BUILDER":
                    raise ValueError(
                        "GraphBuilder nodes cannot be run remotely. Please set run_remotely=False."
                    )
                # this is a PythonJob
                identifier, _ = build_PythonJob_task(identifier)
            return super().new(identifier, name, uuid, **kwargs)
        if isinstance(identifier, str) and identifier.upper() == "PYTHONJOB":
            # copy the inputs and outputs from the function node to the PythonJob node
            identifier, _ = build_PythonJob_task(kwargs.pop("function"))
            return super().new(identifier, name, uuid, **kwargs)
        if isinstance(identifier, str) and identifier.upper() == "SHELLJOB":
            # copy the inputs and outputs from the function node to the SHELLJob node
            identifier, _, links = build_ShellJob_task(
                nodes=kwargs.get("nodes", {}),
                outputs=kwargs.get("outputs", None),
                parser_outputs=kwargs.pop("parser_outputs", None),
            )
            node = super().new(identifier, name, uuid, **kwargs)
            # make links between the tasks
            node.set(links)
            return node
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
