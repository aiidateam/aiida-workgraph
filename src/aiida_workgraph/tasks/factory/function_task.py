from copy import deepcopy
from aiida_workgraph.task import Task
from node_graph.utils import list_to_dict
from aiida_workgraph.orm.mapping import type_mapping
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from node_graph.executor import NodeExecutor
from aiida_workgraph.utils import validate_task_inout


class DecoratedFunctionTaskFactory:
    """
    A factory to create specialized subclasses of Task,
    embedding the 'ndata' (i.e., all relevant data).
    """

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
                "node_class": {
                    "module_path": "aiida_workgraph.task",
                    "callable_name": "Task",
                },
                "group_inputs": group_inputs or [],
                "group_outputs": group_outputs or [],
            },
            "properties": properties,
            "inputs": task_inputs,
            "outputs": task_outputs,
        }
        tdata["executor"] = NodeExecutor.from_callable(func).to_dict()
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        func.identifier = identifier
        return TaskCls

    def __new__(cls, ndata: dict):
        """
        Instead of returning an instance of DecoratedFunctionTaskFactory,
        we return the new Task subclass that uses `ndata`.
        """

        class _DecoratedFunctionTask(Task):
            """A specialized Task with the embedded ndata."""

            _ndata = deepcopy(ndata)

            def __init__(self, *args, **kwargs):
                self.identifier = self._ndata["identifier"]
                self.node_type = self._ndata.get("metadata", {}).get(
                    "node_type", "NORMAL"
                )
                self.catalog = self._ndata.get("metadata", {}).get("catalog", "Others")
                self._error_handlers = self._ndata.get("error_handlers", [])
                super().__init__(*args, **kwargs)
                self.group_inputs = ndata["metadata"].get("group_inputs", [])
                self.group_outputs = ndata["metadata"].get("group_outputs", [])

            def create_properties(self):
                properties = list_to_dict(self._ndata.get("properties") or {})
                for prop in properties.values():
                    self.add_property(
                        prop.pop("identifier", type_mapping["default"]), **prop
                    )

            def create_sockets(self):
                inputs = list_to_dict(self._ndata.get("inputs") or {})
                for inp in inputs.values():
                    if isinstance(inp, str):
                        inp = {"identifier": type_mapping["default"], "name": inp}
                    kwargs = {}
                    if "property_data" in inp:
                        kwargs["property_data"] = inp.pop("property_data")
                    if "sockets" in inp:
                        kwargs["sockets"] = inp.pop("sockets")
                    self.add_input(
                        inp.get("identifier", type_mapping["default"]),
                        name=inp["name"],
                        metadata=inp.get("metadata", {}),
                        link_limit=inp.get("link_limit", 1),
                        **kwargs,
                    )

                outputs = list_to_dict(self._ndata.get("outputs") or {})
                for out in outputs.values():
                    if isinstance(out, str):
                        out = {"identifier": type_mapping["default"], "name": out}
                    self.add_output(
                        out.get("identifier", type_mapping["default"]),
                        name=out["name"],
                        metadata=out.get("metadata", {}),
                    )

            def get_executor(self):
                return self._ndata.get("executor", None)

        return _DecoratedFunctionTask
