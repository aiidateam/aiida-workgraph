from __future__ import annotations
from aiida_workgraph.task import Task
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from aiida_workgraph.orm.mapping import type_mapping
from node_graph.executor import NodeExecutor
import inspect
from aiida.engine.processes.ports import PortNamespace
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from .base import BaseTaskFactory
from aiida_workgraph.utils import validate_task_inout
from aiida_workgraph.tasks.aiida import (
    CalcFunctionTask,
    WorkFunctionTask,
    CalcJobTask,
    WorkChainTask,
)
from node_graph.utils import list_to_dict


task_class_mapping = {
    "CALCFUNCTION": CalcFunctionTask,
    "WORKFUNCTION": WorkFunctionTask,
    "CALCJOB": CalcJobTask,
    "WORKCHAIN": WorkChainTask,
}


def add_input_recursive(
    inputs: List[List[Union[str, Dict[str, Any]]]],
    port: PortNamespace,
    prefix: Optional[str] = None,
    required: bool = True,
) -> List[List[Union[str, Dict[str, Any]]]]:
    """Add input recursively."""
    if prefix is None:
        port_name = port.name
    else:
        port_name = f"{prefix}.{port.name}"
    required = port.required and required
    input_names = [input["name"] for input in inputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port.dynamic:
            link_limit = 1e6
        else:
            link_limit = 1
        if port_name not in input_names:
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": port_name,
                    "link_limit": link_limit,
                    "metadata": {
                        "arg_type": "kwargs",
                        "required": required,
                        "dynamic": port.dynamic,
                    },
                }
            )
        for value in port.values():
            add_input_recursive(inputs, value, prefix=port_name, required=required)
    else:
        if port_name not in input_names:
            # port.valid_type can be a single type or a tuple of types,
            # we only support single type for now
            if isinstance(port.valid_type, tuple) and len(port.valid_type) > 1:
                socket_type = "workgraph.any"
            if isinstance(port.valid_type, tuple) and len(port.valid_type) == 1:
                socket_type = type_mapping.get(port.valid_type[0], "workgraph.any")
            else:
                socket_type = type_mapping.get(port.valid_type, "workgraph.any")
            inputs.append(
                {
                    "identifier": socket_type,
                    "name": port_name,
                    "metadata": {"arg_type": "kwargs", "required": required},
                }
            )
    return inputs


def add_output_recursive(
    outputs: List[List[Union[str, Dict[str, Any]]]],
    port: PortNamespace,
    prefix: Optional[str] = None,
    required: bool = True,
) -> List[List[Union[str, Dict[str, Any]]]]:
    """Add output recursively."""
    if prefix is None:
        port_name = port.name
    else:
        port_name = f"{prefix}.{port.name}"
    required = port.required and required
    output_names = [output["name"] for output in outputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port_name not in output_names:
            outputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": port_name,
                    "metadata": {"required": required, "dynamic": port.dynamic},
                }
            )
        for value in port.values():
            add_output_recursive(outputs, value, prefix=port_name, required=required)
    else:
        if port_name not in output_names:
            outputs.append(
                {
                    "identifier": "workgraph.any",
                    "name": port_name,
                    "metadata": {"required": required},
                }
            )
    return outputs


def get_task_data_from_aiida_component(
    callable: Callable,
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Tuple["Task", Dict[str, Any]]:
    """Register a task from a AiiDA component.
    For example: CalcJob, WorkChain, CalcFunction, WorkFunction."""
    from aiida_workgraph.utils import inspect_aiida_component_type

    tdata = {"metadata": {}}
    inputs = [] if inputs is None else inputs
    outputs = [] if outputs is None else outputs
    task_type = inspect_aiida_component_type(callable)
    if not task_type:
        raise ValueError(f"The callable {callable} is not a valid AiiDA component.")
    spec = callable.spec()
    for _, port in spec.inputs.ports.items():
        add_input_recursive(inputs, port, required=port.required)
    for _, port in spec.outputs.ports.items():
        add_output_recursive(outputs, port, required=port.required)
    # Only check this for calcfunction and workfunction
    if inspect.isfunction(callable) and spec.inputs.dynamic:
        if hasattr(callable.process_class, "_varargs"):
            name = callable.process_class._varargs
        else:
            name = (
                callable.process_class._var_keyword
                or callable.process_class._var_positional
            )
        # if user already defined the var_args in the inputs, skip it
        if name not in [input["name"] for input in inputs]:
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": name,
                    "link_limit": 1e6,
                    "metadata": {"arg_type": "var_kwargs", "dynamic": True},
                }
            )

    tdata["identifier"] = tdata.pop("identifier", callable.__name__)
    tdata["executor"] = NodeExecutor.from_callable(callable).to_dict()
    if task_type.upper() in ["CALCFUNCTION", "AUTO_CALCFUNCTION", "WORKFUNCTION"]:
        outputs = (
            [{"identifier": "workgraph.any", "name": "result"}]
            if not outputs
            else outputs
        )
        tdata["default_name"] = callable.__name__
    # add built-in sockets
    for output in builtin_outputs:
        outputs.append(output.copy())
    for input_data in builtin_inputs:
        inputs.append(input_data.copy())
    final_inputs = {
        "name": "inputs",
        "identifier": "workgraph.namespace",
        "sockets": list_to_dict(inputs),
        "metadata": {"dynamic": spec.inputs.dynamic},
    }
    final_outputs = {
        "name": "outputs",
        "identifier": "workgraph.namespace",
        "sockets": list_to_dict(outputs),
        "metadata": {"dynamic": spec.outputs.dynamic},
    }
    tdata["metadata"]["node_type"] = task_type
    tdata["metadata"]["node_class"] = task_class_mapping.get(task_type.upper(), Task)
    tdata["inputs"] = final_inputs
    tdata["outputs"] = final_outputs
    return tdata


class AiiDAComponentTaskFactory(BaseTaskFactory):
    """Factory for creating tasks from AiiDA components."""

    @classmethod
    def from_aiida_component(
        cls,
        callable: Callable,
        identifier: Optional[str] = None,
        inputs: Optional[List[Union[str, dict]]] = None,
        outputs: Optional[List[Union[str, dict]]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
        additional_data: Optional[Dict[str, Any]] = None,
    ):
        """
        Build the _AiiDAComponentTask subclass from the function
        and the various decorator arguments.
        """

        identifier = identifier or callable.__name__
        error_handlers = error_handlers or []
        inputs = validate_task_inout(inputs or [], "inputs")
        outputs = validate_task_inout(outputs or [], "outputs")

        tdata = get_task_data_from_aiida_component(
            callable, inputs=inputs, outputs=outputs
        )
        tdata["metadata"]["catalog"] = catalog
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        callable.identifier = identifier
        return TaskCls
