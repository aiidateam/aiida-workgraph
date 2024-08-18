from __future__ import annotations
from aiida_workgraph.orm.serializer import serialize_to_aiida_nodes
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.utils import instantiate_process, prepare_inputs
from aiida.manage import manager
from aiida.engine import run_get_node
from aiida.common import InvalidOperation
from aiida.common.log import AIIDA_LOGGER
from aiida.engine.processes import Process, ProcessBuilder
from aiida.orm import ProcessNode
import typing as t
import time

TYPE_SUBMIT_PROCESS = t.Union[Process, t.Type[Process], ProcessBuilder]
LOGGER = AIIDA_LOGGER.getChild("engine.launch")


def prepare_for_workgraph_task(task: dict, kwargs: dict) -> tuple:
    """Prepare the inputs for WorkGraph task"""
    from aiida_workgraph.utils import merge_properties
    from aiida.orm.utils.serialize import deserialize_unsafe

    print("Task type: workgraph.")
    wgdata = deserialize_unsafe(task["executor"]["wgdata"])
    wgdata["name"] = task["name"]
    wgdata["metadata"]["group_outputs"] = task["metadata"]["group_outputs"]
    # update the workgraph data by kwargs
    for task_name, data in kwargs.items():
        # because kwargs is updated using update_nested_dict_with_special_keys
        # which means the data is grouped by the task name
        for socket_name, value in data.items():
            wgdata["tasks"][task_name]["properties"][socket_name]["value"] = value
    # merge the properties
    merge_properties(wgdata)
    metadata = {"call_link_label": task["name"]}
    inputs = {"wg": wgdata, "metadata": metadata}
    return inputs, wgdata


def prepare_for_python_task(task: dict, kwargs: dict, var_kwargs: dict) -> dict:
    """Prepare the inputs for PythonJob"""
    from aiida_workgraph.utils import get_or_create_code
    import os

    print("Task  type: Python.")
    # get the names kwargs for the PythonJob, which are the inputs before _wait
    function_kwargs = {}
    for input in task["inputs"]:
        if input["name"] == "_wait":
            break
        function_kwargs[input["name"]] = kwargs.pop(input["name"], None)
    # if the var_kwargs is not None, we need to pop the var_kwargs from the kwargs
    # then update the function_kwargs if var_kwargs is not None
    if task["metadata"]["var_kwargs"] is not None:
        function_kwargs.pop(task["metadata"]["var_kwargs"], None)
        if var_kwargs:
            # var_kwargs can be AttributeDict if it get data from the previous task output
            if isinstance(var_kwargs, (dict, AttributeDict)):
                function_kwargs.update(var_kwargs)
            # otherwise, it should be a Data node
            elif isinstance(var_kwargs, orm.Data):
                function_kwargs.update(var_kwargs.value)
            else:
                raise ValueError(f"Invalid var_kwargs type: {type(var_kwargs)}")
    # setup code
    code = kwargs.pop("code", None)
    computer = kwargs.pop("computer", None)
    code_label = kwargs.pop("code_label", None)
    code_path = kwargs.pop("code_path", None)
    prepend_text = kwargs.pop("prepend_text", None)
    upload_files = kwargs.pop("upload_files", {})
    new_upload_files = {}
    # change the string in the upload files to SingleFileData, or FolderData
    for key, source in upload_files.items():
        # only alphanumeric and underscores are allowed in the key
        # replace all "." with "_dot_"
        new_key = key.replace(".", "_dot_")
        if isinstance(source, str):
            if os.path.isfile(source):
                new_upload_files[new_key] = orm.SinglefileData(file=source)
            elif os.path.isdir(source):
                new_upload_files[new_key] = orm.FolderData(tree=source)
        elif isinstance(source, (orm.SinglefileData, orm.FolderData)):
            new_upload_files[new_key] = source
        else:
            raise ValueError(f"Invalid upload file type: {type(source)}, {source}")
    #
    if code is None:
        code = get_or_create_code(
            computer=computer if computer else "localhost",
            code_label=code_label if code_label else "python3",
            code_path=code_path if code_path else None,
            prepend_text=prepend_text if prepend_text else None,
        )
    metadata = kwargs.pop("metadata", {})
    metadata.update({"call_link_label": task["name"]})
    # get the source code of the function
    function_name = task["executor"]["function_name"]
    function_source_code = (
        task["executor"]["import_statements"]
        + "\n"
        + task["executor"]["function_source_code_without_decorator"]
    )
    # outputs
    function_outputs = task["outputs"]
    # serialize the kwargs into AiiDA Data
    function_kwargs = serialize_to_aiida_nodes(function_kwargs)
    # transfer the args to kwargs
    inputs = {
        "process_label": f"PythonJob<{task['name']}>",
        "function_source_code": orm.Str(function_source_code),
        "function_name": orm.Str(function_name),
        "code": code,
        "function_kwargs": function_kwargs,
        "upload_files": new_upload_files,
        "function_outputs": orm.List(function_outputs),
        "metadata": metadata,
        **kwargs,
    }
    return inputs


def prepare_for_shell_task(task: dict, kwargs: dict) -> dict:
    """Prepare the inputs for ShellJob"""
    from aiida_shell.launch import prepare_code, convert_nodes_single_file_data
    from aiida.common import lang
    from aiida.orm import AbstractCode

    print("Task  type: ShellJob.")
    command = kwargs.pop("command", None)
    resolve_command = kwargs.pop("resolve_command", False)
    metadata = kwargs.pop("metadata", {})
    # setup code
    if isinstance(command, str):
        computer = (metadata or {}).get("options", {}).pop("computer", None)
        code = prepare_code(command, computer, resolve_command)
    else:
        lang.type_check(command, AbstractCode)
        code = command
    # update the tasks with links
    nodes = convert_nodes_single_file_data(kwargs.pop("nodes", {}))
    # find all keys in kwargs start with "nodes."
    for key in list(kwargs.keys()):
        if key.startswith("nodes."):
            nodes[key[6:]] = kwargs.pop(key)
    metadata.update({"call_link_label": task["name"]})
    inputs = {
        "code": code,
        "nodes": nodes,
        "filenames": kwargs.pop("filenames", {}),
        "arguments": kwargs.pop("arguments", []),
        "outputs": kwargs.pop("outputs", []),
        "parser": kwargs.pop("parser", None),
        "metadata": metadata or {},
    }
    return inputs


# modified from aiida.engine.submit
# do not check the scope of the process
def submit(
    process: TYPE_SUBMIT_PROCESS,
    inputs: dict[str, t.Any] | None = None,
    *,
    wait: bool = False,
    wait_interval: int = 5,
    **kwargs: t.Any,
) -> ProcessNode:
    """Submit the process with the supplied inputs to the daemon immediately returning control to the interpreter.

    .. warning: this should not be used within another process. Instead, there one should use the ``submit`` method of
        the wrapping process itself, i.e. use ``self.submit``.

    .. warning: submission of processes requires ``store_provenance=True``.

    :param process: the process class, instance or builder to submit
    :param inputs: the input dictionary to be passed to the process
    :param wait: when set to ``True``, the submission will be blocking and wait for the process to complete at which
        point the function returns the calculation node.
    :param wait_interval: the number of seconds to wait between checking the state of the process when ``wait=True``.
    :param kwargs: inputs to be passed to the process. This is an alternative to the positional ``inputs`` argument.
    :return: the calculation node of the process
    """
    inputs = prepare_inputs(inputs, **kwargs)

    runner = manager.get_manager().get_runner()
    assert runner.persister is not None, "runner does not have a persister"
    assert runner.controller is not None, "runner does not have a controller"

    process_inited = instantiate_process(runner, process, **inputs)

    # If a dry run is requested, simply forward to `run`, because it is not compatible with `submit`. We choose for this
    # instead of raising, because in this way the user does not have to change the launcher when testing. The same goes
    # for if `remote_folder` is present in the inputs, which means we are importing an already completed calculation.
    if process_inited.metadata.get("dry_run", False) or "remote_folder" in inputs:
        _, node = run_get_node(process_inited)
        return node

    if not process_inited.metadata.store_provenance:
        raise InvalidOperation("cannot submit a process with `store_provenance=False`")

    runner.persister.save_checkpoint(process_inited)
    process_inited.close()

    # Do not wait for the future's result, because in the case of a single worker this would cock-block itself
    runner.controller.continue_process(process_inited.pid, nowait=False, no_reply=True)
    node = process_inited.node

    if not wait:
        return node

    while not node.is_terminated:
        LOGGER.report(
            f"Process<{node.pk}> has not yet terminated, current state is `{node.process_state}`. "
            f"Waiting for {wait_interval} seconds."
        )
        time.sleep(wait_interval)

    return node
