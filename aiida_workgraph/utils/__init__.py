from typing import Any, Dict, Optional, Union, Callable
from aiida.engine.processes import Process
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.engine.runners import Runner


def get_executor(data: Dict[str, Any]) -> Union[Process, Any]:
    """Import executor from path and return the executor and type."""
    import importlib
    from aiida.plugins import CalculationFactory, WorkflowFactory, DataFactory

    is_pickle = data.get("is_pickle", False)
    type = data.get("type", "function")
    if is_pickle:
        import cloudpickle as pickle

        try:
            executor = pickle.loads(data["executor"])
        except Exception as e:
            print("Error in loading executor: ", e)
            executor = None
    else:
        if type == "WorkflowFactory":
            executor = WorkflowFactory(data["name"])
        elif type == "CalculationFactory":
            executor = CalculationFactory(data["name"])
        elif type == "DataFactory":
            executor = DataFactory(data["name"])
        else:
            module = importlib.import_module("{}".format(data.get("path", "")))
            executor = getattr(module, data["name"])

    return executor, type


def create_data_node(executor: orm.Data, args: list, kwargs: dict) -> orm.Node:
    """Create an AiiDA data node from the executor and args and kwargs."""
    from aiida import orm

    # print("Create data node: ", executor, args, kwargs)
    extras = kwargs.pop("extras", {})
    repository_metadata = kwargs.pop("repository_metadata", {})
    if issubclass(executor, (orm.BaseType, orm.Dict)):
        data_node = executor(*args)
    else:
        data_node = executor(*args)
        data_node.base.attributes.set_many(kwargs)
    data_node.base.extras.set_many(extras)
    data_node.base.repository.repository_metadata = repository_metadata
    data_node.store()
    return data_node


def get_nested_dict(d: Dict, name: str) -> Any:
    """
    name = "base.pw.parameters"
    """
    keys = name.split(".")
    current = d
    for key in keys:
        if key not in current:
            raise ValueError(f"Context variable {name} not found.")
        current = current[key]
    return current


def update_nested_dict(d: Dict[str, Any], key: str, value: Any) -> None:
    """
    d = {}
    key = "base.pw.parameters"
    value = 2
    will give:
    d = {"base": {"pw": {"parameters": 2}}
    """
    keys = key.split(".")
    current = d
    current = {} if current is None else current
    for k in keys[:-1]:
        current = current.setdefault(k, {})
    current[keys[-1]] = value


def is_empty(value: Any) -> bool:
    """Check if the provided value is an empty collection."""
    import numpy as np

    if isinstance(value, np.ndarray):
        return value.size == 0
    elif isinstance(value, (dict, list, set, tuple)):
        return not value
    return False


def update_nested_dict_with_special_keys(d: Dict[str, Any]) -> Dict[str, Any]:
    """Remove None and empty value"""
    d = {k: v for k, v in d.items() if v is not None and not is_empty(v)}
    #
    special_keys = [k for k in d.keys() if "." in k]
    for key in special_keys:
        value = d.pop(key)
        update_nested_dict(d, key, value)
    return d


def merge_properties(wgdata: Dict[str, Any]) -> None:
    """Merge sub properties to the root properties.
    {
        "base.pw.parameters": 2,
        "base.pw.code": 1,
    }
    after merge:
    {"base": {"pw": {"parameters": 2,
                    "code": 1}}
    So that no "." in the key name.
    """
    for name, task in wgdata["tasks"].items():
        for key, prop in task["properties"].items():
            if "." in key and prop["value"] not in [None, {}]:
                root, key = key.split(".", 1)
                update_nested_dict(
                    task["properties"][root]["value"], key, prop["value"]
                )
                prop["value"] = None


def generate_node_graph(pk: int) -> Any:
    from aiida.tools.visualization import Graph
    from aiida import orm

    graph = Graph()
    calc_node = orm.load_node(pk)
    graph.recurse_ancestors(calc_node, annotate_links="both")
    graph.recurse_descendants(calc_node, annotate_links="both")
    return graph.graphviz


def build_task_link(wgdata: Dict[str, Any]) -> None:
    """Create links for tasks.
    Create the links for task inputs using:
    1) workgraph links
    2) if it is a graph builder graph, expose the group inputs and outputs
    sockets.
    """
    # reset task input links
    for name, task in wgdata["tasks"].items():
        for input in task["inputs"]:
            input["links"] = []
        for output in task["outputs"]:
            output["links"] = []
    for link in wgdata["links"]:
        to_socket = [
            socket
            for socket in wgdata["tasks"][link["to_node"]]["inputs"]
            if socket["name"] == link["to_socket"]
        ][0]
        from_socket = [
            socket
            for socket in wgdata["tasks"][link["from_node"]]["outputs"]
            if socket["name"] == link["from_socket"]
        ][0]
        to_socket["links"].append(link)
        from_socket["links"].append(link)


def get_dict_from_builder(builder: Any) -> Dict:
    """Transform builder to pure dict."""
    from aiida.engine.processes.builder import ProcessBuilderNamespace

    if isinstance(builder, ProcessBuilderNamespace):
        return {k: get_dict_from_builder(v) for k, v in builder.items()}
    else:
        return builder


def get_workgraph_data(process: Union[int, orm.Node]) -> Optional[Dict[str, Any]]:
    """Get the workgraph data from the process node."""
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import load_node

    if isinstance(process, int):
        process = load_node(process)
    wgdata = process.base.extras.get("_workgraph", None)
    if wgdata is None:
        return
    wgdata = deserialize_unsafe(wgdata)
    return wgdata


def get_parent_workgraphs(pk: int) -> list:
    """Get the list of parent workgraphs.
    Use aiida incoming links to find the parent workgraphs.
    the parent workgraph is the workgraph that has a link (type CALL_WORK) to the current workgraph.
    """
    from aiida import orm
    from aiida.common.links import LinkType

    node = orm.load_node(pk)
    parent_workgraphs = [[node.process_label, node.pk]]
    links = node.base.links.get_incoming(link_type=LinkType.CALL_WORK).all()
    print(links)
    if len(links) > 0:
        parent_workgraphs.extend(get_parent_workgraphs(links[0].node.pk))
    print(parent_workgraphs)
    return parent_workgraphs


def get_processes_latest(pk: int) -> Dict[str, Dict[str, Union[int, str]]]:
    """Get the latest info of all tasks from the process."""
    import aiida
    from aiida.orm.utils.serialize import deserialize_unsafe

    process = aiida.orm.load_node(pk)
    tasks = {}
    for key in process.base.extras.keys():
        if key.startswith("_task_state"):
            name = key[12:]
            state = deserialize_unsafe(process.base.extras.get(key))
            task_process = deserialize_unsafe(
                process.base.extras.get(f"_task_process_{name}")
            )
            if task_process:
                tasks[name] = {
                    "pk": task_process.pk,
                    "state": state,
                    "ctime": task_process.ctime,
                    "mtime": task_process.mtime,
                }
            else:
                tasks[name] = {
                    "pk": None,
                    "state": state,
                    "ctime": None,
                    "mtime": None,
                }
    return tasks


def get_or_create_code(
    computer: str = "localhost",
    code_label: str = "python3",
    code_path: str = None,
    prepend_text: str = "",
):
    """Try to load code, create if not exit."""
    from aiida.orm.nodes.data.code.installed import InstalledCode

    try:
        return orm.load_code(f"{code_label}@{computer}")
    except NotExistent:
        description = f"Python code on computer: {computer}"
        computer = orm.load_computer(computer)
        code_path = code_path or code_label
        code = InstalledCode(
            computer=computer,
            label=code_label,
            description=description,
            filepath_executable=code_path,
            default_calc_job_plugin="workgraph.python",
            prepend_text=prepend_text,
        )

        code.store()
        return code


def serialize_pythonjob_properties(wgdata):
    """Serialize the PythonJob properties."""
    from aiida_workgraph.orm.serializer import general_serializer

    for _, task in wgdata["tasks"].items():
        if not task["metadata"]["node_type"].upper() == "PYTHONJOB":
            continue
        # get the names kwargs for the PythonJob, which are the inputs before _wait
        input_kwargs = []
        for input in task["inputs"]:
            if input["name"] == "_wait":
                break
            input_kwargs.append(input["name"])
        for name in input_kwargs:
            prop = task["properties"][name]
            # if value is not None, not {}
            if not (
                prop["value"] is None
                or isinstance(prop["value"], dict)
                and prop["value"] == {}
            ):
                prop["value"] = general_serializer(prop["value"])


def generate_bash_to_create_python_env(
    name: str,
    pip: list = None,
    conda: dict = None,
    modules: list = None,
    python_version: str = None,
    variables: dict = None,
    shell: str = "posix",
):
    """
    Generates a bash script for creating or updating a Python environment on a remote computer.
    If python_version is None, it uses the Python version from the local environment.
    Conda is a dictionary that can include 'channels' and 'dependencies'.
    """
    import sys

    pip = pip or []
    conda_channels = conda.get("channels", []) if conda else []
    conda_dependencies = conda.get("dependencies", []) if conda else []
    # Determine the Python version from the local environment if not provided
    local_python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    desired_python_version = (
        python_version if python_version is not None else local_python_version
    )

    # Start of the script
    script = "#!/bin/bash\n\n"

    # Load modules if provided
    if modules:
        script += "# Load specified system modules\n"
        for module in modules:
            script += f"module load {module}\n"

    # Conda shell hook initialization for proper conda activation
    script += "# Initialize Conda for this shell\n"
    script += f'eval "$(conda shell.{shell} hook)"\n'

    script += "# Setup the Python environment\n"
    script += "if ! conda info --envs | grep -q ^{name}$; then\n"
    script += "    # Environment does not exist, create it\n"
    if conda_dependencies:
        dependencies_string = " ".join(conda_dependencies)
        script += f"    conda create -y -n {name} python={desired_python_version} {dependencies_string}\n"
    else:
        script += f"    conda create -y -n {name} python={desired_python_version}\n"
    script += "fi\n"
    if conda_channels:
        script += "EXISTING_CHANNELS=$(conda config --show channels)\n"
        script += "for CHANNEL in " + " ".join(conda_channels) + ";\n"
        script += "do\n"
        script += '    if ! echo "$EXISTING_CHANNELS" | grep -q $CHANNEL; then\n'
        script += "        conda config --prepend channels $CHANNEL\n"
        script += "    fi\n"
        script += "done\n"
    script += f"conda activate {name}\n"

    # Install pip packages
    if pip:
        script += f"pip install {' '.join(pip)}\n"

    # Set environment variables
    if variables:
        for var, value in variables.items():
            script += f"export {var}='{value}'\n"

    # End of the script
    script += "echo 'Environment setup is complete.'\n"

    return script


def create_conda_env(
    computer: Union[str, orm.Computer],
    name: str,
    pip: list = None,
    conda: list = None,
    modules: list = None,
    python_version: str = None,
    variables: dict = None,
    shell: str = "posix",
) -> tuple:
    """Test that there is no unexpected output from the connection."""
    # Execute a command that should not return any error, except ``NotImplementedError``
    # since not all transport plugins implement remote command execution.
    from aiida.common.exceptions import NotExistent
    from aiida import orm

    user = orm.User.collection.get_default()
    if isinstance(computer, str):
        computer = orm.load_computer(computer)
    try:
        authinfo = computer.get_authinfo(user)
    except NotExistent:
        raise f"Computer<{computer.label}> is not yet configured for user<{user.email}>"

    scheduler = authinfo.computer.get_scheduler()
    transport = authinfo.get_transport()

    script = generate_bash_to_create_python_env(
        name, pip, conda, modules, python_version, variables, shell
    )
    with transport:
        scheduler.set_transport(transport)
        try:
            retval, stdout, stderr = transport.exec_command_wait(script)
        except NotImplementedError:
            return (
                True,
                f"Skipped, remote command execution is not implemented for the "
                f"`{computer.transport_type}` transport plugin",
            )

        if retval != 0:
            return (
                False,
                f"The command `echo -n` returned a non-zero return code ({retval})",
            )

        template = """
We detected an error while creating the environemnt on the remote computer, as shown between the bars
=============================================================================================
{}
=============================================================================================
Please check!
    """
        if stderr:
            return False, template.format(stderr)

        if stdout:
            # the last line is the echo 'Environment setup is complete.'
            if not stdout.strip().endswith("Environment setup is complete."):
                return False, template.format(stdout)
            else:
                return True, "Environment setup is complete."

    return True, None


def create_and_pause_process(
    runner: Runner = None,
    process_class: Callable = None,
    inputs: dict = None,
    state_msg: str = "",
) -> Process:
    from aiida.engine.utils import instantiate_process

    process_inited = instantiate_process(runner, process_class, **inputs)
    process_inited.pause(msg=state_msg)
    process_inited.runner.persister.save_checkpoint(process_inited)
    process_inited.close()
    runner.controller.continue_process(process_inited.pid, nowait=True, no_reply=True)
    return process_inited


def recursive_to_dict(attr_dict):
    """
    Recursively convert an AttributeDict to a standard dictionary.
    """
    from aiida.common import AttributeDict
    from plumpy.utils import AttributesFrozendict

    if isinstance(attr_dict, (AttributesFrozendict, AttributeDict)):
        return {k: recursive_to_dict(v) for k, v in attr_dict.items()}
    else:
        return attr_dict
