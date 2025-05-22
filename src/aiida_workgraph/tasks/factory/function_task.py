from aiida_workgraph.orm.mapping import type_mapping
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from node_graph.executor import NodeExecutor
from aiida_workgraph.utils import validate_task_inout
from .base import BaseTaskFactory
from node_graph.utils import list_to_dict


class DecoratedFunctionTaskFactory(BaseTaskFactory):
    """A factory to create specialized subclasses of Task from functions."""

    default_task_type = "Normal"

    @classmethod
    def from_function(
        cls,
        func: Callable,
        identifier: Optional[str] = None,
        task_type: str = None,
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Union[str, dict]]] = None,
        outputs: Optional[List[Union[str, dict]]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
        additional_data: Optional[Dict[str, Any]] = None,
        node_class: Optional[Callable] = None,
    ):
        """
        Build the _DecoratedFunctionTask subclass from the function
        and the various decorator arguments.
        """
        from node_graph.decorator import generate_input_sockets

        node_class = node_class or cls.default_base_class

        task_type = task_type or cls.default_task_type
        identifier = identifier or func.__name__
        inputs = inputs or []
        properties = properties or []
        # at least one output is required
        task_outputs = outputs or [{"identifier": "workgraph.any", "name": "result"}]
        error_handlers = error_handlers or []
        inputs = validate_task_inout(inputs, "inputs")
        task_outputs = validate_task_inout(task_outputs, "outputs")
        task_inputs = generate_input_sockets(
            func, inputs, properties, type_mapping=type_mapping
        )
        # Mark function inputs and outputs
        task_outputs = {
            "name": "outputs",
            "identifier": node_class.SocketPool.any,
            "sockets": list_to_dict(task_outputs),
        }
        for out in task_outputs["sockets"].values():
            out.setdefault("metadata", {})
            out["metadata"]["function_socket"] = True
        # add built-in sockets
        for input_data in builtin_inputs:
            task_inputs["sockets"][input_data["name"]] = input_data.copy()
        for output in builtin_outputs:
            task_outputs["sockets"][output["name"]] = output.copy()

        tdata = {
            "identifier": identifier,
            "metadata": {
                "node_type": task_type,
                "catalog": catalog,
            },
            "properties": properties,
            "inputs": task_inputs,
            "outputs": task_outputs,
            "error_handlers": error_handlers,
        }
        tdata["executor"] = NodeExecutor.from_callable(func).to_dict()
        if node_class:
            tdata["metadata"]["node_class"] = node_class
        tdata["default_name"] = func.__name__
        additional_data = additional_data or {}
        tdata.update(additional_data)

        NodeCls = cls(tdata)
        return NodeCls
