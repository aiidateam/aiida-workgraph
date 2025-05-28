from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from aiida_workgraph.orm.mapping import type_mapping
from node_graph.executor import NodeExecutor


def generate_tdata(
    func: Callable,
    identifier: str,
    inputs: List[Tuple[str, str]],
    outputs: List[Tuple[str, str]],
    properties: List[Tuple[str, str]],
    catalog: str,
    task_type: str,
    additional_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate task data for creating a task."""
    from node_graph.decorator import generate_input_sockets

    task_inputs = generate_input_sockets(
        func, inputs, properties, type_mapping=type_mapping
    )
    for input_data in task_inputs:
        input_data.setdefault("metadata", {})
        input_data["metadata"]["function_socket"] = True
    task_outputs = outputs
    # add built-in sockets
    for output in builtin_outputs:
        task_outputs.append(output.copy())
    for input_data in builtin_inputs:
        task_inputs.append(input_data.copy())
    tdata = {
        "identifier": identifier,
        "metadata": {
            "task_type": task_type,
            "catalog": catalog,
            "node_class": {
                "module_path": "aiida_workgraph.task",
                "callable_name": "Task",
            },
        },
        "properties": properties,
        "inputs": task_inputs,
        "outputs": task_outputs,
    }
    tdata["executor"] = NodeExecutor.from_callable(func).to_dict()
    if additional_data:
        tdata.update(additional_data)
    return tdata
