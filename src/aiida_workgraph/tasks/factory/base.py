from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple, Union
from node_graph.utils import list_to_dict
from aiida_workgraph.orm.mapping import type_mapping
from aiida_workgraph.task import Task
from node_graph.executor import NodeExecutor


class BaseTaskFactory:
    """
    A base factory to create specialized subclasses of Task,
    embedding the 'ndata' (i.e., all relevant data).
    """

    @classmethod
    def create_task(
        cls,
        identifier: str,
        task_type: str,
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Union[str, dict]]] = None,
        outputs: Optional[List[Union[str, dict]]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
        additional_data: Optional[Dict[str, Any]] = None,
        executor: Optional[Dict[str, Any]] = None,
    ):
        """
        Create the _Task subclass from provided task details.
        """
        properties = properties or []
        inputs = inputs or []
        outputs = outputs or []
        error_handlers = error_handlers or []

        tdata = {
            "identifier": identifier,
            "metadata": {
                "node_type": task_type,
                "catalog": catalog,
            },
            "properties": properties,
            "inputs": inputs,
            "outputs": outputs,
            "error_handlers": error_handlers,
            "executor": executor or {},
        }
        additional_data = additional_data or {}
        tdata.update(additional_data)

        return cls._create_task_class(tdata)

    def __new__(cls, ndata: dict):
        class _TaskFactory(Task):
            """A specialized Task with the embedded ndata."""

            _ndata = deepcopy(ndata)

            def __init__(self, *args, **kwargs):
                self.identifier = self._ndata["identifier"]
                self.node_type = self._ndata.get("metadata", {}).get(
                    "node_type", "NORMAL"
                )
                self.catalog = self._ndata.get("metadata", {}).get("catalog", "Others")
                self._executor = NodeExecutor(**self._ndata["executor"])
                self._error_handlers = self._ndata.get("error_handlers", [])
                super().__init__(*args, **kwargs)
                self.group_inputs = ndata["metadata"].get("group_inputs", [])
                self.group_outputs = ndata["metadata"].get("group_outputs", [])

            def create_properties(self):
                properties = list_to_dict(self._ndata.get("properties", {}))
                for prop in properties.values():
                    self.add_property(
                        prop.pop("identifier", type_mapping["default"]), **prop
                    )

            def create_sockets(self):
                inputs = list_to_dict(self._ndata.get("inputs", {}))
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

                outputs = list_to_dict(self._ndata.get("outputs", {}))
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

        return _TaskFactory
