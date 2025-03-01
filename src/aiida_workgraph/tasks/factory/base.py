from node_graph.utils import list_to_dict
from aiida_workgraph.orm.mapping import type_mapping
from aiida_workgraph.task import Task
import importlib


class BaseTaskFactory:
    """
    A base factory to create specialized subclasses of Task,
    embedding the 'ndata' (i.e., all relevant data).
    """

    def __new__(cls, ndata: dict):

        ndata.setdefault("metadata", {})
        BaseClass = ndata["metadata"].get("node_class", Task)
        if isinstance(BaseClass, dict):
            module_path = BaseClass["module_path"]
            callable_name = BaseClass["callable_name"]
            module = importlib.import_module(module_path)
            BaseClass = getattr(module, callable_name)

        class _TaskFactory(BaseClass):
            """A specialized Task with the embedded ndata."""

            _ndata = ndata
            is_dynamic: bool = True

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
                properties = list_to_dict(self._ndata.get("properties", {}))
                for prop in properties.values():
                    prop.setdefault("identifier", type_mapping["default"])
                    self.add_property(**prop)

            def create_sockets(self):
                inputs = list_to_dict(self._ndata.get("inputs", {}))
                for inp in inputs.values():
                    if isinstance(inp, str):
                        inp = {"identifier": type_mapping["default"], "name": inp}
                    kwargs = {}
                    if "property_data" in inp:
                        kwargs["property_data"] = inp.get("property_data", {})
                    if "sockets" in inp:
                        kwargs["sockets"] = inp.get("sockets", None)
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

            def get_metadata(self):
                metadata = super().get_metadata()
                metadata["node_class"] = {
                    "module_path": BaseClass.__module__,
                    "callable_name": BaseClass.__name__,
                }
                metadata["factory_class"] = {
                    "module_path": BaseTaskFactory.__module__,
                    "callable_name": BaseTaskFactory.__name__,
                }
                return metadata

        return _TaskFactory
