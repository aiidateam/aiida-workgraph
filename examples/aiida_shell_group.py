from aiida_worktree import node, WorkTree, build_node
from aiida import load_profile
import pathlib
import shlex
import typing as t

from aiida.common import lang
from aiida.orm import AbstractCode, Data, ProcessNode
from aiida.parsers import Parser

from aiida_shell.launch import prepare_code, convert_nodes_single_file_data

load_profile()

shelljob = build_node({"path": "aiida_shell.calculations.shell.ShellJob"})


@node.group(
    outputs=[
        ["shelljob1.stdout", "stdout"],
        ["shelljob1.stderr", "stderr"],
        ["shelljob1.retrieved", "retrieved"],
        ["shelljob1.remote_folder", "remote_folder"],
    ]
)
def launch_shell_job(  # noqa: PLR0913
    command: str | AbstractCode,
    arguments: list[str] | str | None = None,
    nodes: t.Mapping[str, str | pathlib.Path | Data] | None = None,
    filenames: dict[str, str] | None = None,
    outputs: list[str] | None = None,
    parser: t.Callable[[Parser, pathlib.Path], dict[str, Data]] | str | None = None,
    metadata: dict[str, t.Any] | None = None,
) -> tuple[dict[str, Data], ProcessNode]:
    """Launch a :class:`aiida_shell.ShellJob` job for the given command.

    :param command: The shell command to run. Should be the relative command name, e.g., ``date``. An ``AbstractCode``
        instance will be automatically created for this command if it doesn't already exist. Alternatively, a pre-
        configured ``AbstractCode`` instance can be passed directly.
    :param arguments: Optional list of command line arguments optionally containing placeholders for input nodes. The
        arguments can also be specified as a single string. In this case, it will be split into separate parameters
        using ``shlex.split``.
    :param nodes: A dictionary of ``Data`` nodes whose content is to replace placeholders in the ``arguments`` list.
    :param filenames: Optional dictionary of explicit filenames to use for the ``nodes`` to be written to ``dirpath``.
    :param outputs: Optional list of relative filenames that should be captured as outputs.
    :param parser: Optional callable that can implement custom parsing logic of produced output files. Alternatively,
        a complete entry point, i.e. a string of the form ``{entry_point_group}:{entry_point_name}`` pointing to such a
        callable.
    :param metadata: Optional dictionary of metadata inputs to be passed to the ``ShellJob``.
    :param submit: Boolean, if ``True`` will submit the job to the daemon instead of running in current interpreter.
    :raises TypeError: If the value specified for ``metadata.options.computer`` is not a ``Computer``.
    :raises ValueError: If the absolute path of the command on the computer could not be determined.
    :returns: The tuple of results dictionary and ``ProcessNode``, or just the ``ProcessNode`` if ``submit=True``. The
        results dictionary intentionally doesn't include the ``retrieved`` and ``remote_folder`` outputs as they are
        generated for each ``CalcJob`` and typically are not of interest to a user running ``launch_shell_job``. In
        order to not confuse them, these nodes are omitted, but they can always be accessed through the node.
    """
    computer = (metadata or {}).get("options", {}).pop("computer", None)

    if isinstance(command, str):
        code = prepare_code(command, computer)
    else:
        lang.type_check(command, AbstractCode)
        code = command

    if isinstance(arguments, str):
        arguments = shlex.split(arguments)
    else:
        lang.type_check(arguments, list, allow_none=True)

    inputs = {
        "code": code,
        "nodes": convert_nodes_single_file_data(nodes or {}),
        "filenames": filenames,
        "arguments": arguments,
        "outputs": outputs,
        "parser": parser,
        "metadata": metadata or {},
    }

    wt = WorkTree(name="test_aiida_shell")
    shelljob1 = wt.nodes.new(shelljob, "shelljob1")
    shelljob1.set(inputs)

    return wt


@node()
def generate_nodes(data):
    """Prepare the nodes"""
    return {"pdb": data}


# Create a worktree
wt = WorkTree(name="test_aiida_shell")
job1 = wt.nodes.new(launch_shell_job, command="pdb_fetch", arguments=["1brs"])
job2 = wt.nodes.new(launch_shell_job, command="pdb_selchain", arguments=["-A,D {pdb}"])
job3 = wt.nodes.new(launch_shell_job, command="pdb_delhetatm", arguments=["{pdb}"])
job4 = wt.nodes.new(launch_shell_job, command="pdb_tidy", arguments=["{pdb}"])
generate_nodes1 = wt.nodes.new(generate_nodes)
generate_nodes2 = wt.nodes.new(generate_nodes)
generate_nodes3 = wt.nodes.new(generate_nodes)
wt.links.new(job1.outputs["stdout"], generate_nodes1.inputs[0])
wt.links.new(generate_nodes1.outputs[0], job2.inputs["nodes"])
wt.links.new(job2.outputs["stdout"], generate_nodes2.inputs[0])
wt.links.new(generate_nodes2.outputs[0], job3.inputs["nodes"])
wt.links.new(job3.outputs["stdout"], generate_nodes3.inputs[0])
wt.links.new(generate_nodes3.outputs[0], job4.inputs["nodes"])
wt.submit()
