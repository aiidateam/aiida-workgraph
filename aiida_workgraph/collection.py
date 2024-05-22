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
            build_node_from_callable,
            build_PythonJob_node,
        )

        # build the node on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_node_from_callable(identifier)
            if kwargs.pop("run_remotely", False):
                # this is a PythonJob
                identifier = build_PythonJob_node(identifier)
        if isinstance(identifier, str) and identifier.upper() == "PYTHONJOB":
            # copy the inputs and outputs from the function node to the PythonJob node
            identifier = build_PythonJob_node(kwargs.pop("function"))
        # Call the original new method
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
