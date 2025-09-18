from __future__ import annotations
from node_graph.node_spec import NodeSpec
from aiida_workgraph.task import SpecTask


class SubGraphTask(SpecTask):
    """Task created from WorkGraph."""

    identifier = "workgraph.workgraph_task"
    name = "SubGraphTask"
    node_type = "Normal"
    catalog = "Builtins"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._subgraph = None

    @property
    def subgraph(self):
        from aiida_workgraph import WorkGraph
        from copy import deepcopy

        if not self._subgraph:
            graph_data = deepcopy(self.get_executor().graph_data)
            self._subgraph = WorkGraph.from_dict(graph_data)
        return self._subgraph

    @property
    def tasks(self):
        return self.subgraph.tasks

    @property
    def links(self):
        return self.subgraph.links

    def prepare_for_subgraph_task(self, kwargs: dict) -> tuple:
        """Prepare the inputs for SubGraph task"""
        # update the subgraph inputs by the kwargs
        for name, data in kwargs.items():
            input_socket = self.subgraph.inputs[name]
            input_socket._set_socket_value(data)
        # merge the properties
        metadata = {"call_link_label": self.name}
        inputs = self.subgraph.to_engine_inputs(metadata=metadata)
        return inputs

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        inputs = self.prepare_for_subgraph_task(kwargs)

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


def _build_subgraph_task_nodespec(
    graph: "WorkGraph",
    name: str | None = None,
) -> NodeSpec:
    from node_graph.executor import SafeExecutor

    meta = {
        "node_type": "SubGraph",
    }

    return NodeSpec(
        identifier=name or graph.name,
        inputs=graph._inputs,
        outputs=graph._outputs,
        executor=SafeExecutor.from_graph(graph),
        base_class=SubGraphTask,
        metadata=meta,
    )
