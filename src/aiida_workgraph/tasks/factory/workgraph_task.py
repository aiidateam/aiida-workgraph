from typing import TYPE_CHECKING
from .base import BaseTaskFactory
from aiida_workgraph.config import builtin_inputs, builtin_outputs

if TYPE_CHECKING:
    from aiida_workgraph import WorkGraph


class WorkGraphTaskFactory(BaseTaskFactory):
    """A factory to create Task from WorkGraph."""

    @classmethod
    def create_task(
        cls,
        workgraph: "WorkGraph",
    ):
        tdata = {"metadata": {"node_type": "workgraph"}}
        inputs = []
        outputs = []
        group_outputs = []
        # add all the inputs/outputs from the tasks in the workgraph
        builtin_input_names = [input["name"] for input in builtin_inputs]
        builtin_output_names = [output["name"] for output in builtin_outputs]

        for task in workgraph.tasks:
            # inputs
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": f"{task.name}",
                }
            )
            for socket in task.inputs:
                if socket._name in builtin_input_names:
                    continue
                inputs.append(
                    {
                        "identifier": socket._identifier,
                        "name": f"{task.name}.{socket._name}",
                    }
                )
            # outputs
            outputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": f"{task.name}",
                }
            )
            for socket in task.outputs:
                if socket._name in builtin_output_names:
                    continue
                outputs.append(
                    {
                        "identifier": socket._identifier,
                        "name": f"{task.name}.{socket._name}",
                    }
                )
                group_outputs.append(
                    {
                        "name": f"{task.name}.{socket._name}",
                        "from": f"{task.name}.{socket._name}",
                    }
                )
        # add built-in sockets
        for output in builtin_outputs:
            outputs.append(output.copy())
        for input in builtin_inputs:
            inputs.append(input.copy())
        tdata["metadata"]["node_class"] = {
            "module_path": "aiida_workgraph.tasks.builtins",
            "callable_name": "WorkGraphTask",
        }
        tdata["inputs"] = inputs
        tdata["outputs"] = outputs
        tdata["identifier"] = workgraph.name
        # get graph_data from the workgraph
        graph_data = workgraph.prepare_inputs()["workgraph_data"]
        executor = {
            "module_path": "aiida_workgraph.engine.workgraph",
            "callable_name": "WorkGraphEngine",
            "graph_data": graph_data,
        }
        tdata["metadata"]["group_outputs"] = group_outputs
        tdata["executor"] = executor

        TaskCls = cls(tdata)
        return TaskCls
