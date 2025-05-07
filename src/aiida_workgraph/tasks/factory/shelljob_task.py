from typing import List, Optional, Union, Dict, Any
from .base import BaseTaskFactory
from aiida_shell import ShellJob
from .aiida_task import AiiDAComponentTaskFactory
from aiida_workgraph.utils import validate_task_inout
from aiida_workgraph import Task
from node_graph.executor import NodeExecutor

additional_inputs = [
    {"identifier": "workgraph.any", "name": "command"},
    {"identifier": "workgraph.any", "name": "resolve_command"},
]
additional_outputs = [
    {"identifier": "workgraph.any", "name": "stdout"},
    {"identifier": "workgraph.any", "name": "stderr"},
]


def prepare_for_shell_task(inputs: dict) -> dict:
    """Prepare the inputs for ShellJob"""
    from aiida_shell.launch import prepare_shell_job_inputs
    import inspect

    # Retrieve the signature of `prepare_shell_job_inputs` to determine expected input parameters.
    signature = inspect.signature(prepare_shell_job_inputs)
    aiida_shell_input_keys = signature.parameters.keys()

    # Iterate over all WorkGraph `inputs`, and extract the ones which are expected by `prepare_shell_job_inputs`
    inputs_aiida_shell_subset = {
        key: inputs[key] for key in inputs.keys() if key in aiida_shell_input_keys
    }
    # if parser in inputs, and the parser is a dict
    parser = inputs_aiida_shell_subset.get("parser", None)
    if isinstance(parser, dict):
        inputs_aiida_shell_subset["parser"] = NodeExecutor(**parser).executor

    try:
        aiida_shell_inputs = prepare_shell_job_inputs(**inputs_aiida_shell_subset)
    except ValueError:
        raise

    # We need to remove the original input-keys, as they might be offending for the call to `launch_shell_job`
    # E.g., `inputs` originally can contain `command`, which gets, however, transformed to #
    # `code` by `prepare_shell_job_inputs`
    for key in inputs_aiida_shell_subset.keys():
        inputs.pop(key)

    # Finally, we update the original `inputs` with the modified ones from the call to `prepare_shell_job_inputs`
    inputs = {**inputs, **aiida_shell_inputs}

    inputs.setdefault("metadata", {})
    return inputs


class ShellJobTask(Task):
    """Task with AiiDA calcfunction/workfunction as executor."""

    identifier = "workgraph.shelljob"
    name = "shelljob"
    node_type = "ShellJob"
    catalog = "AIIDA"

    def serialize_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Overwrite the serialize_data method to handle the parser function."""
        import inspect

        prop = data["inputs"]["sockets"]["parser"]["property"]
        if prop["value"] is not None:
            if inspect.isfunction(prop["value"]):
                prop["value"] = NodeExecutor.from_callable(prop["value"]).to_dict()

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process

        inputs = prepare_for_shell_task(kwargs)
        inputs["metadata"].update({"call_link_label": self.name})
        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                ShellJob,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(ShellJob, **inputs)
            state = "RUNNING"
        process.label = self.name

        return process, state


class ShellJobTaskFactory(BaseTaskFactory):
    """A factory to create ShellJobTask."""

    @classmethod
    def create_class(cls, inputs: dict) -> Union[Task, None]:
        return ShellJobTaskFactory.create_task(
            outputs=inputs.get("outputs", None),
            parser_outputs=inputs.pop("parser_outputs", None),
        )

    @classmethod
    def create_task(
        cls,
        outputs: Optional[List[Union[str, dict]]] = None,
        parser_outputs: Optional[List[Union[str, dict]]] = None,
    ):
        from aiida_shell.parsers.shell import ShellParser

        outputs = outputs or []
        parser_outputs = parser_outputs or []
        parser_outputs = validate_task_inout(parser_outputs, "parser_outputs")
        TaskCls = AiiDAComponentTaskFactory.from_aiida_component(ShellJob)
        tdata = TaskCls._ndata
        for output in additional_outputs:
            tdata["outputs"]["sockets"][output["name"]] = output.copy()

        outputs = [
            {
                "identifier": "workgraph.any",
                "name": ShellParser.format_link_label(output),
            }
            for output in outputs
        ]
        outputs.extend(parser_outputs)
        for output in outputs:
            if output["name"] not in tdata["outputs"]["sockets"]:
                tdata["outputs"]["sockets"][output["name"]] = output.copy()
        #
        tdata["identifier"] = "ShellJob"
        for input_data in additional_inputs:
            tdata["inputs"]["sockets"][input_data["name"]] = input_data.copy()
        tdata["metadata"]["node_type"] = "SHELLJOB"
        tdata["metadata"]["node_class"] = ShellJobTask

        TaskCls = cls(tdata)
        return TaskCls
