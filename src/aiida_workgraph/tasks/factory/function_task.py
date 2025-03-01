from aiida_workgraph.orm.mapping import type_mapping
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from node_graph.executor import NodeExecutor
from aiida_workgraph.utils import validate_task_inout
from .base import BaseTaskFactory


class DecoratedFunctionTaskFactory(BaseTaskFactory):
    """A factory to create specialized subclasses of Task from functions."""

    @classmethod
    def from_function(
        cls,
        func: Callable,
        identifier: Optional[str] = None,
        task_type: str = "Normal",
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Union[str, dict]]] = None,
        outputs: Optional[List[Union[str, dict]]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
        group_inputs: List[Tuple[str, str]] = None,
        group_outputs: List[Tuple[str, str]] = None,
        additional_data: Optional[Dict[str, Any]] = None,
        node_class: Optional[Callable] = None,
    ):
        """
        Build the _DecoratedFunctionTask subclass from the function
        and the various decorator arguments.
        """
        from node_graph.decorator import generate_input_sockets

        identifier = identifier or func.__name__
        inputs = inputs or []
        properties = properties or []
        task_outputs = outputs or []
        error_handlers = error_handlers or []
        inputs = validate_task_inout(inputs, "inputs")
        task_outputs = validate_task_inout(task_outputs, "outputs")
        task_inputs = generate_input_sockets(
            func, inputs, properties, type_mapping=type_mapping
        )
        # Mark function inputs and outputs
        for input in task_inputs:
            input.setdefault("metadata", {})
            input["metadata"]["is_function_input"] = True
        for out in task_outputs:
            out.setdefault("metadata", {})
            out["metadata"]["is_function_output"] = True
        # add built-in sockets
        for input in builtin_inputs:
            task_inputs.append(input.copy())
        for output in builtin_outputs:
            task_outputs.append(output.copy())
        tdata = {
            "identifier": identifier,
            "metadata": {
                "node_type": task_type,
                "catalog": catalog,
                "group_inputs": group_inputs or [],
                "group_outputs": group_outputs or [],
            },
            "properties": properties,
            "inputs": task_inputs,
            "outputs": task_outputs,
            "error_handlers": error_handlers,
        }
        tdata["executor"] = NodeExecutor.from_callable(func).to_dict()
        if node_class:
            tdata["metadata"]["node_class"] = node_class
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        func.identifier = identifier
        return TaskCls
