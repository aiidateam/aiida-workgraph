from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple, Union
from node_graph.utils import list_to_dict
from aiida_workgraph.orm.mapping import type_mapping
from aiida_workgraph.task import Task


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

    @staticmethod
    def _create_task_class(ndata: dict):
        """Dynamically create a specialized Task class embedding the given data."""

        class _Task(Task):
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

            def create_properties(self):
                properties = list_to_dict(self._ndata.get("properties") or {})
                for prop in properties.values():
                    self.add_property(
                        prop.pop("identifier", type_mapping["default"]), **prop
                    )

            def create_sockets(self):
                for key in ["inputs", "outputs"]:
                    items = list_to_dict(self._ndata.get(key) or {})
                    for item in items.values():
                        if isinstance(item, str):
                            item = {"identifier": type_mapping["default"], "name": item}
                        kwargs = {"metadata": item.get("metadata", {})}
                        if key == "inputs":
                            kwargs["link_limit"] = item.get("link_limit", 1)
                            self.add_input(
                                item.get("identifier", type_mapping["default"]),
                                name=item["name"],
                                **kwargs
                            )
                        else:
                            self.add_output(
                                item.get("identifier", type_mapping["default"]),
                                name=item["name"],
                                **kwargs
                            )

            def get_executor(self):
                return self._ndata.get("executor", None)

        return _Task
