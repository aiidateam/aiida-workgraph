from typing import Any, Dict, Optional, Union, Callable, List
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


def get_nested_dict(d: Dict, name: str, allow_none: bool = False) -> Any:
    """
    name = "base.pw.parameters"
    """
    keys = name.split(".")
    current = d
    for key in keys:
        if key not in current:
            if allow_none:
                return None
            else:
                raise ValueError(f"{name} not exist in {d}")
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


def get_pythonjob_data(tdata: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process the task data dictionary for a PythonJob.
    It load the orignal Python data from the AiiDA Data node for the
    args and kwargs of the function.

    Args:
        tdata (Dict[str, Any]): The input data dictionary.

    Returns:
        Dict[str, Any]: The processed data dictionary.
    """
    for name in tdata["metadata"]["args"]:
        if tdata["properties"][name]["value"] is None:
            continue
        if name in tdata["properties"]:
            tdata["properties"][name]["value"] = tdata["properties"][name][
                "value"
            ].value
    for name in tdata["metadata"]["kwargs"]:
        # all the kwargs are after computer is the input for the PythonJob, should be AiiDA Data node
        if tdata["properties"][name]["value"] is None:
            continue
        if name == "computer":
            break
        if name in tdata["properties"]:
            tdata["properties"][name]["value"] = tdata["properties"][name][
                "value"
            ].value
    return tdata


def serialize_workgraph_data(wgdata: Dict[str, Any]) -> Dict[str, Any]:
    from aiida.orm.utils.serialize import serialize

    for name, task in wgdata["tasks"].items():
        wgdata["tasks"][name] = serialize(task)
    wgdata["error_handlers"] = serialize(wgdata["error_handlers"])


def get_workgraph_data(process: Union[int, orm.Node]) -> Optional[Dict[str, Any]]:
    """Get the workgraph data from the process node."""
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import load_node

    if isinstance(process, int):
        process = load_node(process)
    wgdata = process.base.extras.get("_workgraph", None)
    if wgdata is None:
        return
    for name, task in wgdata["tasks"].items():
        wgdata["tasks"][name] = deserialize_unsafe(task)
        if wgdata["tasks"][name]["metadata"]["node_type"].upper() == "PYTHONJOB":
            get_pythonjob_data(wgdata["tasks"][name])
    wgdata["error_handlers"] = deserialize_unsafe(wgdata["error_handlers"])
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


def get_processes_latest(
    pk: int, node_name: str = None
) -> Dict[str, Dict[str, Union[int, str]]]:
    """Get the latest info of all tasks from the process."""
    import aiida
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import QueryBuilder
    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    tasks = {}
    node_names = [node_name] if node_name else []
    if node_name:
        projections = [
            f"extras._task_state_{node_name}",
            f"extras._task_process_{node_name}",
        ]
    else:
        projections = []
        process = aiida.orm.load_node(pk)
        node_names = [
            key[12:]
            for key in process.base.extras.keys()
            if key.startswith("_task_state")
        ]
        projections = [f"extras._task_state_{name}" for name in node_names]
        projections.extend([f"extras._task_process_{name}" for name in node_names])
    qb = QueryBuilder()
    qb.append(WorkGraphEngine, filters={"id": pk}, project=projections)
    # print("projections: ", projections)
    results = qb.all()
    # change results to dict
    results = dict(zip(projections, results[0]))
    # print("results: ", results)
    for name in node_names:
        state = results[f"extras._task_state_{name}"]
        task_process = deserialize_unsafe(results[f"extras._task_process_{name}"])
        tasks[name] = {
            "pk": task_process.pk if task_process else None,
            "process_type": task_process.process_type if task_process else "",
            "state": state,
            "ctime": task_process.ctime if task_process else None,
            "mtime": task_process.mtime if task_process else None,
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


def get_required_imports(func):
    """Retrieve type hints and the corresponding module"""
    from typing import get_type_hints, _SpecialForm

    type_hints = get_type_hints(func)
    imports = {}

    def add_imports(type_hint):
        if isinstance(
            type_hint, _SpecialForm
        ):  # Handle special forms like Any, Union, Optional
            module_name = "typing"
            type_name = type_hint._name or str(type_hint)
        elif hasattr(
            type_hint, "__origin__"
        ):  # This checks for higher-order types like List, Dict
            module_name = type_hint.__module__
            type_name = getattr(type_hint, "_name", None) or getattr(
                type_hint.__origin__, "__name__", None
            )
            for arg in getattr(type_hint, "__args__", []):
                if arg is type(None):  # noqa: E721
                    continue
                add_imports(arg)  # Recursively add imports for each argument
        elif hasattr(type_hint, "__module__"):
            module_name = type_hint.__module__
            type_name = type_hint.__name__
        else:
            return  # If no module or origin, we can't import it, e.g., for literals

        if module_name not in imports:
            imports[module_name] = set()
        if type_name is not None:
            imports[module_name].add(type_name)

    for _, type_hint in type_hints.items():
        add_imports(type_hint)

    return imports


def serialize_function(func: Callable) -> Dict[str, Any]:
    """Serialize a function for storage or transmission."""
    import inspect
    import textwrap
    import cloudpickle as pickle

    # we need save the source code explicitly, because in the case of jupyter notebook,
    # the source code is not saved in the pickle file
    source_code = inspect.getsource(func)
    # Split the source into lines for processing
    source_code_lines = source_code.split("\n")
    function_source_code = "\n".join(source_code_lines)
    # Find the first line of the actual function definition
    for i, line in enumerate(source_code_lines):
        if line.strip().startswith("def "):
            break
    function_source_code_without_decorator = "\n".join(source_code_lines[i:])
    function_source_code_without_decorator = textwrap.dedent(
        function_source_code_without_decorator
    )
    # we also need to include the necessary imports for the types used in the type hints.
    try:
        required_imports = get_required_imports(func)
    except Exception as e:
        required_imports = {}
        print(f"Failed to get required imports for function {func.__name__}: {e}")
    # Generate import statements
    import_statements = "\n".join(
        f"from {module} import {', '.join(types)}"
        for module, types in required_imports.items()
    )
    return {
        "executor": pickle.dumps(func),
        "type": "function",
        "is_pickle": True,
        "function_name": func.__name__,
        "function_source_code": function_source_code,
        "function_source_code_without_decorator": function_source_code_without_decorator,
        "import_statements": import_statements,
    }


def workgraph_to_short_json(
    wgdata: Dict[str, Union[str, List, Dict]]
) -> Dict[str, Union[str, Dict]]:
    """Export a workgraph to a rete js editor data."""
    wgdata_short = {
        "name": wgdata["name"],
        "uuid": wgdata["uuid"],
        "state": wgdata["state"],
        "nodes": {},
        "links": wgdata["links"],
    }
    #
    for name, task in wgdata["tasks"].items():
        # Add required inputs to nodes
        inputs = [
            input
            for input in task["inputs"]
            if input["name"] in task["metadata"]["args"]
        ]
        wgdata_short["nodes"][name] = {
            "label": task["name"],
            "inputs": inputs,
            "outputs": [],
            "position": task["position"],
        }
    # Add links to nodes
    for link in wgdata["links"]:
        wgdata_short["nodes"][link["to_node"]]["inputs"].append(
            {
                "name": link["to_socket"],
            }
        )
        wgdata_short["nodes"][link["from_node"]]["outputs"].append(
            {
                "name": link["from_socket"],
            }
        )
    return wgdata_short
