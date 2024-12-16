from aiida_workgraph.task import Task
from aiida_workgraph.utils.inspect_aiida_components import (
    get_task_data_from_aiida_component,
)
from aiida_workgraph.orm.mapping import type_mapping
from copy import deepcopy
from node_graph.utils import list_to_dict
from node_graph.executor import NodeExecutor


class AiiDAProcessTask(Task):
    """Task with AiiDA process as executor."""

    identifier = "workgraph.aiida_process"
    name = "aiida_process"
    node_type = "Process"
    catalog = "AIIDA"

    def __init__(self, executor: NodeExecutor, **kwargs):

        task_data = self.inspect_executor(executor)
        self.task_data = task_data
        super().__init__(executor=executor, **kwargs)
        self.node_type = self.task_data["metadata"]["task_type"]

    def inspect_executor(self, executor) -> dict:
        metadata = executor.get("metadata", {})
        inputs = metadata.pop("inputs", [])
        outputs = metadata.pop("outputs", [])
        task_data = get_task_data_from_aiida_component(
            tdata=deepcopy(executor),
            inputs=deepcopy(inputs),
            outputs=deepcopy(outputs),
        )
        return task_data

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()

        inputs = list_to_dict(self.task_data.get("inputs", {}))
        for input in inputs.values():
            if isinstance(input, str):
                input = {"identifier": type_mapping["default"], "name": input}
            kwargs = {}
            if "property_data" in input:
                kwargs["property_data"] = input.pop("property_data")
            self.add_input(
                input.get("identifier", type_mapping["default"]),
                name=input["name"],
                metadata=input.get("metadata", {}),
                link_limit=input.get("link_limit", 1),
                **kwargs,
            )
        outputs = list_to_dict(self.task_data.get("outputs", {}))
        for output in outputs.values():
            if isinstance(output, str):
                output = {"identifier": type_mapping["default"], "name": output}
            identifier = output.get("identifier", type_mapping["default"])
            self.add_output(
                identifier, name=output["name"], metadata=output.get("metadata", {})
            )
