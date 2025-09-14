from __future__ import annotations

from dataclasses import replace
from typing import Any, Dict, List, Optional, Union, Callable
import inspect
from aiida_shell import ShellJob
from aiida_shell.launch import prepare_shell_job_inputs
from node_graph.node_spec import NodeSpec
from node_graph.executor import RuntimeExecutor
from node_graph.socket_spec import SocketSpec, merge_specs
from aiida_workgraph.socket_spec import from_aiida_process, namespace
from aiida_workgraph.task import SpecTask, TaskHandle
from aiida import orm


class ShellJobTask(SpecTask):
    """Runtime for ShellJob nodes.

    This class is referenced by NodeSpec.base_class_path so the engine can import
    it and call `execute`.
    """

    identifier = "workgraph.shelljob"
    name = "shelljob"
    node_type = "SHELLJOB"
    catalog = "AIIDA"

    def serialize_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Overwrite the serialize_data method to handle the parser function."""
        import inspect

        parser = data["inputs"].get("parser")
        if parser is not None:
            if inspect.isfunction(parser):
                data["inputs"]["parser"] = RuntimeExecutor.from_callable(
                    parser
                ).to_dict()

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        """Submit/launch the AiiDA ShellJob.

        - Translates friendly inputs (command, resolve_command, parser, ...)
          using `prepare_shell_job_inputs`.
        - Submits or runs under the engine's runner.
        """
        from aiida_workgraph.utils import create_and_pause_process

        kwargs = dict(kwargs or {})

        # Detect and translate aiida-shell convenience arguments
        signature = inspect.signature(prepare_shell_job_inputs)
        aiida_shell_keys = signature.parameters.keys()

        subset = {k: kwargs[k] for k in list(kwargs) if k in aiida_shell_keys}

        parser = subset.get("parser", None)
        if isinstance(parser, dict) and {"module_path", "callable_name"} <= set(parser):
            # already a Executor dict -> build executor instance
            subset["parser"] = RuntimeExecutor(**parser).callable
        elif inspect.isfunction(parser):
            subset["parser"] = parser

        if subset:
            if "command" in subset:
                subset["command"] = (
                    subset["command"].value
                    if isinstance(subset["command"], orm.Str)
                    else subset["command"]
                )
            if "resolve_command" in subset:
                subset["resolve_command"] = (
                    subset["resolve_command"].value
                    if isinstance(subset["resolve_command"], orm.Bool)
                    else subset["resolve_command"]
                )
            if "arguments" in subset:
                subset["arguments"] = (
                    subset["arguments"].get_list()
                    if isinstance(subset["arguments"], orm.List)
                    else subset["arguments"]
                )
            prepared = prepare_shell_job_inputs(**subset)
            # drop original keys so they won't clash with launch kwargs
            for k in subset.keys():
                kwargs.pop(k, None)
            # merge translated inputs
            kwargs.update(prepared)

        # metadata
        md = kwargs.setdefault("metadata", {})
        md.setdefault("call_link_label", self.name)

        if getattr(self, "action", None) == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                ShellJob,
                kwargs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(ShellJob, **kwargs)
            state = "RUNNING"
        process.label = self.name
        return process, state


def _build_shelljob_nodespec(
    *,
    identifier: Optional[str] = None,
    outputs: Optional[SocketSpec | List[str]] = None,
    parser_outputs: Optional[SocketSpec | List[str]] = None,
) -> NodeSpec:
    """Create a `NodeSpec` for a ShellJob, augmenting inputs/outputs as needed.

    - Start from AiiDA Process spec inference
    - Add inputs: command, resolve_command
    - Ensure stdout/stderr outputs exist
    - Optionally add user-declared outputs and parser_outputs (as leaf-any)
    """
    from aiida_workgraph.socket_spec import validate_socket_data
    from aiida_shell.parsers.shell import ShellParser

    outputs = validate_socket_data(outputs)
    parser_outputs = validate_socket_data(parser_outputs)

    in_spec, out_spec = from_aiida_process(ShellJob)

    # Add additional inputs
    additions_in = namespace(command=Any, resolve_command=bool)
    in_spec = merge_specs(in_spec, additions_in)

    # Ensure stdout/stderr outputs
    additions_out = namespace(stdout=Any, stderr=Any)
    out_spec = merge_specs(out_spec, additions_out)

    # add extra outputs requested by user
    if outputs:
        # make sure the key are AiiDA compatible
        fields = {
            ShellParser.format_link_label(key): value
            for key, value in outputs.fields.items()
        }
        outputs = replace(outputs, fields=fields)
        out_spec = merge_specs(out_spec, outputs)

    if parser_outputs:
        out_spec = merge_specs(out_spec, parser_outputs)

    exec_payload = RuntimeExecutor.from_callable(ShellJob)

    return NodeSpec(
        identifier=identifier or "ShellJob",
        catalog="AIIDA",
        inputs=in_spec,
        outputs=out_spec,
        executor=exec_payload,
        base_class=ShellJobTask,
        metadata={"node_type": "SHELLJOB"},
    )


# Public factory used by users inside a WorkGraph


def shelljob(
    *,
    command: str,
    arguments: Optional[List[str]] = None,
    nodes: Optional[Dict[str, Any]] = None,
    filenames: Optional[Dict[str, str]] = None,
    outputs: Optional[List[Union[str, Dict[str, Any]]]] = None,
    parser: Optional[Callable] = None,
    parser_outputs: Optional[SocketSpec | List[str]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    resolve_command: bool = True,
):
    """Create a ShellJob node in the active WorkGraph and return its outputs handle.

    Usage:
        with WorkGraph(name="test_shell_date_with_arguments") as wg:
            outs = shelljob(command="date", arguments=["--iso-8601"])  # returns handle
            wg.run()
    """
    spec = _build_shelljob_nodespec(outputs=outputs, parser_outputs=parser_outputs)

    handle = TaskHandle(spec)
    return handle(
        command=command,
        arguments=arguments,
        nodes=nodes,
        filenames=filenames,
        outputs=outputs,
        parser=parser,
        metadata=metadata,
        resolve_command=resolve_command,
    )
