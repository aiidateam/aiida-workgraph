from typing import TYPE_CHECKING
from .base import BaseTaskFactory
from aiida_workgraph.config import builtin_inputs, builtin_outputs
from aiida_workgraph import Task

if TYPE_CHECKING:
    from aiida_workgraph import WorkGraph


class WorkGraphTask(Task):
    """Task created from WorkGraph."""

    identifier = "workgraph.workgraph_task"
    name = "WorkGraphTask"
    node_type = "Normal"
    catalog = "Builtins"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._workgraph = None

    @property
    def workgraph(self):
        from aiida_workgraph import WorkGraph
        from copy import deepcopy

        if not self._workgraph:
            graph_data = deepcopy(self.get_executor()["graph_data"])
            self._workgraph = WorkGraph.from_dict(graph_data)
        return self._workgraph

    @property
    def tasks(self):
        return self.workgraph.tasks

    @property
    def links(self):
        return self.workgraph.links

    def prepare_for_workgraph_task(self, kwargs: dict) -> tuple:
        """Prepare the inputs for WorkGraph task"""

        wgdata = self.get_executor()["graph_data"]
        wgdata["name"] = self.name
        wgdata["metadata"]["group_outputs"] = self.metadata["group_outputs"]
        # update the workgraph data by kwargs
        for task_name, data in kwargs.items():
            # because kwargs is updated using update_nested_dict_with_special_keys
            # which means the data is grouped by the task name
            for socket_name, value in data.items():
                input = wgdata["tasks"][task_name]["inputs"]["sockets"][socket_name]
                if input["identifier"] == "workgraph.namespace":
                    input["value"] = value
                else:
                    input["property"]["value"] = value
        # merge the properties
        metadata = {"call_link_label": self.name}
        inputs = {"workgraph_data": wgdata, "metadata": metadata}
        return inputs

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        inputs = self.prepare_for_workgraph_task(kwargs)

        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                WorkGraphEngine,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(WorkGraphEngine, **inputs)
            state = "RUNNING"
        process.label = self.name

        return process, state


class WorkGraphTaskFactory(BaseTaskFactory):
    """A factory to create Task from WorkGraph."""

    @classmethod
    def create_task(
        cls,
        workgraph: "WorkGraph",
    ):
        tdata = {"metadata": {"node_type": "workgraph"}}
        inputs = {"name": "inputs", "identifier": "workgraph.namespace", "sockets": {}}
        outputs = {
            "name": "outputs",
            "identifier": "workgraph.namespace",
            "sockets": {},
        }
        group_outputs = []
        # add all the inputs/outputs from the tasks in the workgraph
        # builtin_input_names = [input["name"] for input in builtin_inputs]
        builtin_output_names = [output["name"] for output in builtin_outputs]

        for task in workgraph.tasks:
            # inputs
            data = task.inputs._to_dict()
            data["name"] = task.name
            inputs["sockets"][task.name] = data
            # outputs
            data = task.outputs._to_dict()
            data["name"] = task.name
            outputs["sockets"][task.name] = data
            for socket in task.outputs:
                if socket._name in builtin_output_names:
                    continue
                group_outputs.append(
                    {
                        "name": f"{task.name}.{socket._name}",
                        "from": f"{task.name}.{socket._name}",
                    }
                )
        # add built-in sockets
        for input in builtin_inputs:
            inputs["sockets"][input["name"]] = input.copy()
        for output in builtin_outputs:
            outputs["sockets"][output["name"]] = output.copy()
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
        tdata["metadata"]["node_class"] = WorkGraphTask
        tdata["executor"] = executor

        TaskCls = cls(tdata)
        return TaskCls
