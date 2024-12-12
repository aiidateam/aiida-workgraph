from __future__ import annotations

from typing import Any, Dict, List, Optional, Union, Tuple
import inspect
from aiida.engine.processes.ports import PortNamespace
from aiida_workgraph.utils import build_callable
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from aiida_workgraph.orm.mapping import type_mapping


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
        if port_name not in input_names:
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": port_name,
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
    tdata: Dict[str, Any],
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Tuple["Task", Dict[str, Any]]:
    """Register a task from a AiiDA component.
    For example: CalcJob, WorkChain, CalcFunction, WorkFunction."""
    import importlib
    from aiida_workgraph.utils import inspect_aiida_component_type

    tdata.setdefault("metadata", {})
    inputs = [] if inputs is None else inputs
    outputs = [] if outputs is None else outputs
    if "module_path" in tdata:
        module = importlib.import_module("{}".format(tdata.get("module_path", "")))
        executor = getattr(module, tdata["callable_name"])
    else:
        executor = tdata["callable"]
    task_type = inspect_aiida_component_type(executor)
    if not task_type:
        raise ValueError(f"The callable {executor} is not a valid AiiDA component.")
    spec = executor.spec()
    for _, port in spec.inputs.ports.items():
        add_input_recursive(inputs, port, required=port.required)
    for _, port in spec.outputs.ports.items():
        add_output_recursive(outputs, port, required=port.required)
    # Only check this for calcfunction and workfunction
    if inspect.isfunction(executor) and spec.inputs.dynamic:
        if hasattr(executor.process_class, "_varargs"):
            name = executor.process_class._varargs
        else:
            name = (
                executor.process_class._var_keyword
                or executor.process_class._var_positional
            )
        # if user already defined the var_args in the inputs, skip it
        if name not in [input["name"] for input in inputs]:
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": name,
                    "metadata": {"arg_type": "var_kwargs", "dynamic": True},
                }
            )

    tdata["identifier"] = tdata.pop("identifier", executor.__name__)
    tdata["executor"] = build_callable(executor)
    if task_type.upper() in ["CALCFUNCTION", "WORKFUNCTION"]:
        outputs = (
            [{"identifier": "workgraph.any", "name": "result"}]
            if not outputs
            else outputs
        )
    # add built-in sockets
    for output in builtin_outputs:
        outputs.append(output.copy())
    for input in builtin_inputs:
        inputs.append(input.copy())
    tdata["metadata"]["node_class"] = {
        "module_path": "aiida_workgraph.task",
        "callable_name": "Task",
    }
    tdata["metadata"]["task_type"] = task_type
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    return tdata
