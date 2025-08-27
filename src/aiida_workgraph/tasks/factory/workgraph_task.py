from aiida_workgraph import Task


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
        # update the workgraph inputs by the kwargs
        for name, data in kwargs.items():
            input_socket = self.workgraph.inputs[name]
            input_socket._set_socket_value(data)
        # merge the properties
        metadata = {"call_link_label": self.name}
        inputs = self.workgraph.prepare_inputs(metadata=metadata)
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
