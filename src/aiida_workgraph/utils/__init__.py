from __future__ import annotations

from typing import Any, Dict, Optional, Union, Callable, List
from aiida.engine.processes import Process
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.engine.runners import Runner
from aiida_workgraph.config import task_types
from aiida.engine import CalcJob, WorkChain
from aiida_pythonjob import PythonJob
from aiida_shell.calculations.shell import ShellJob
import inspect


def inspect_aiida_component_type(executor: Callable) -> str:
    task_type = None
    if isinstance(executor, type):
        if executor == PythonJob:
            task_type = "PYTHONJOB"
        elif executor == ShellJob:
            task_type = "SHELLJOB"
        elif issubclass(executor, CalcJob):
            task_type = task_types[CalcJob]
        elif issubclass(executor, WorkChain):
            task_type = task_types[WorkChain]
    elif inspect.isfunction(executor):
        if getattr(executor, "node_class", False):
            task_type = task_types[executor.node_class]
    return task_type


def get_nested_dict(d: Dict, name: str, **kwargs) -> Any:
    """Get the value from a nested dictionary.
    If default is provided, return the default value if the key is not found.
    Otherwise, raise ValueError.
    For example:
    d = {"base": {"pw": {"parameters": 2}}}
    name = "base.pw.parameters"
    """
    from aiida.orm.utils.managers import NodeLinksManager

    keys = name.split(".")
    current = d
    for key in keys:
        if key not in current:
            if "default" in kwargs:
                return kwargs.get("default")
            else:
                if isinstance(current, dict):
                    avaiable_keys = current.keys()
                elif isinstance(current, NodeLinksManager):
                    avaiable_keys = list(current._get_keys())
                else:
                    avaiable_keys = []
                raise ValueError(f"{name} not exist. Available keys: {avaiable_keys}")
        current = current[key]
    return current


def merge_dicts(dict1: Any, dict2: Any) -> Any:
    """Recursively merges two dictionaries."""
    for key, value in dict2.items():
        if key in dict1 and isinstance(dict1[key], dict) and isinstance(value, dict):
            # Recursively merge dictionaries
            dict1[key] = merge_dicts(dict1[key], value)
        else:
            # Overwrite or add the key
            dict1[key] = value
    return dict1


def update_nested_dict(
    base: Optional[Dict[str, Any]], key_path: str, value: Any
) -> None:
    """
    Update or create a nested dictionary structure based on a dotted key path.

    This function allows updating a nested dictionary or creating one if `d` is `None`.
    Given a dictionary and a key path (e.g., "base.pw.parameters"), it will traverse
    or create the necessary nested structure to set the provided value at the specified
    key location. If intermediate dictionaries do not exist, they will be created.
    If the resulting dictionary is empty, it is set to `None`.

    Args:
        base (Dict[str, Any] | None): The dictionary to update, which can be `None`.
                                   If `None`, an empty dictionary will be created.
        key (str): A dotted key path string representing the nested structure.
        value (Any): The value to set at the specified key.

    Example:
        base = None
        key = "scf.pw.parameters"
        value = 2
        After running:
            update_nested_dict(d, key, value)
        The result will be:
            base = {"scf": {"pw": {"parameters": 2}}}

    Edge Case:
        If the resulting dictionary is empty after the update, it will be set to `None`.

    """

    if base is None:
        base = {}
    keys = key_path.split(".")
    current_key = keys[0]
    if len(keys) == 1:
        # Base case: Merge dictionaries or set the value directly.
        if isinstance(base.get(current_key), dict) and isinstance(value, dict):
            base[current_key] = merge_dicts(base[current_key], value)
        else:
            base[current_key] = value
    else:
        # Recursive case: Ensure the key exists and is a dictionary, then recurse.
        if current_key not in base or not isinstance(base[current_key], dict):
            base[current_key] = {}
        base[current_key] = update_nested_dict(
            base[current_key], ".".join(keys[1:]), value
        )

    return base


def update_nested_dict_with_special_keys(data: Dict[str, Any]) -> Dict[str, Any]:
    """Update the nested dictionary with special keys like "base.pw.parameters"."""
    # Remove None

    data = {k: v for k, v in data.items() if v is not None}
    #
    special_keys = [k for k in data.keys() if "." in k]
    for key in special_keys:
        value = data.pop(key)
        update_nested_dict(data, key, value)
    return data


def generate_node_graph(
    pk: int, output: str = None, width: str = "100%", height: str = "600px"
) -> Any:
    """Generate the node graph for the given node pk.
    If in Jupyter, return the graphviz object.
    Otherwise, save the graph to an HTML file.
    """

    from aiida.tools.visualization import Graph
    from aiida import orm
    from IPython.display import IFrame
    import pathlib
    from .svg_to_html import svg_to_html

    in_jupyter = False
    try:
        from IPython import get_ipython

        if get_ipython() is not None:
            in_jupyter = True
    except NameError:
        pass

    graph = Graph()
    calc_node = orm.load_node(pk)
    graph.recurse_ancestors(calc_node, annotate_links="both")
    graph.recurse_descendants(calc_node, annotate_links="both")
    g = graph.graphviz
    if not in_jupyter:
        html_content = svg_to_html(g._repr_image_svg_xml(), width, height)
        if output is None:
            pathlib.Path("html").mkdir(exist_ok=True)
            output = f"html/node_graph_{pk}.html"
        with open(output, "w") as f:
            f.write(html_content)
        return IFrame(output, width=width, height=height)
    return g


def get_dict_from_builder(builder: Any) -> Dict:
    """Transform builder to pure dict."""
    from aiida.engine.processes.builder import ProcessBuilderNamespace

    if isinstance(builder, ProcessBuilderNamespace):
        return {k: get_dict_from_builder(v) for k, v in builder.items()}
    else:
        return builder


def get_workgraph_data(process: Union[int, orm.Node]) -> Optional[Dict[str, Any]]:
    """Get the workgraph data from the process node."""
    from aiida_workgraph.orm.utils import deserialize_safe
    from aiida.orm import load_node

    if isinstance(process, int):
        process = load_node(process)
    wgdata = process.workgraph_data
    task_executors = process.task_executors
    if wgdata is None:
        return
    for name, task in wgdata["tasks"].items():
        wgdata["tasks"][name] = deserialize_safe(task)
        wgdata["tasks"][name]["executor"] = task_executors.get(name)
    wgdata["error_handlers"] = process.workgraph_error_handlers
    wgdata["context"] = deserialize_safe(wgdata["context"])
    return wgdata


def get_parent_workgraphs(pk: int) -> List[List[str, int]]:
    """Get the list of parent workgraphs.
    Use aiida incoming links to find the parent workgraphs.
    the parent workgraph is the workgraph that has a link (type CALL_WORK) to the current workgraph.
    """
    from aiida import orm
    from aiida.common.links import LinkType

    node = orm.load_node(pk)
    parent_workgraphs = [[node.process_label, node.pk]]
    links = node.base.links.get_incoming(link_type=LinkType.CALL_WORK).all()
    if len(links) > 0:
        parent_workgraphs.extend(get_parent_workgraphs(links[0].node.pk))
    return parent_workgraphs


def get_processes_latest(
    pk: int, task_name: str = None, item_type: str = "task"
) -> Dict[str, Dict[str, Union[int, str]]]:
    """Get the latest info of all tasks from the process."""
    import aiida
    from aiida_workgraph.orm.utils import deserialize_safe

    tasks = {}
    if pk is None:
        return tasks
    node = aiida.orm.load_node(pk)
    if item_type == "called_process":
        # fetch the process that called by the workgraph
        for link in node.base.links.get_outgoing().all():
            if isinstance(link.node, aiida.orm.ProcessNode):
                tasks[f"{link.node.process_label}_{link.node.pk}"] = {
                    "pk": link.node.pk,
                    "process_type": link.node.process_type,
                    "state": link.node.process_state.value,
                    "ctime": link.node.ctime,
                    "mtime": link.node.mtime,
                }
    elif item_type == "task":
        task_states = node.task_states
        task_processes = node.task_processes
        task_names = [task_name] if task_name else task_states.keys()
        for name in task_names:
            state = task_states[name]
            task_process = deserialize_safe(task_processes.get(name, ""))
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
        description = f"Code on computer: {computer}"
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


def create_and_pause_process(
    runner: Runner = None,
    process_class: Callable = None,
    inputs: dict = None,
    state_msg: str = "",
) -> Process:
    from aiida.engine.utils import instantiate_process

    process_inited = instantiate_process(runner, process_class, **inputs)
    process_inited.pause(state_msg)
    process_inited.runner.persister.save_checkpoint(process_inited)
    process_inited.close()
    runner.controller.continue_process(process_inited.pid, nowait=True, no_reply=True)
    return process_inited


def get_raw_value(identifier, value: Any) -> Any:
    """Get the raw value from a Data node."""
    if identifier in [
        "workgraph.int",
        "workgraph.float",
        "workgraph.string",
        "workgraph.bool",
        "workgraph.aiida_int",
        "workgraph.aiida_float",
        "workgraph.aiida_string",
        "workgraph.aiida_bool",
    ]:
        if value is not None and isinstance(value, orm.Data):
            return value.value
        else:
            return value
    elif (
        identifier == "workgraph.aiida_structure"
        and value is not None
        and isinstance(value, orm.StructureData)
    ):
        content = value.backend_entity.attributes
        content["node_type"] = value.node_type
        return content
    elif isinstance(value, orm.Data):
        content = value.backend_entity.attributes
        content["node_type"] = value.node_type
        return content


def process_properties(task: Dict) -> Dict:
    """Extract raw values."""
    result = {}
    for name, prop in task["properties"].items():
        identifier = prop["identifier"]
        value = prop.get("value")
        result[name] = {
            "identifier": identifier,
            "value": get_raw_value(identifier, value),
        }
    #
    for name, input in task["inputs"]["sockets"].items():
        if input.get("property"):
            prop = input["property"]
            identifier = prop["identifier"]
            value = prop.get("value")
            result[name] = {
                "identifier": identifier,
                "value": get_raw_value(identifier, value),
            }

    return result


def workgraph_to_short_json(
    wgdata: Dict[str, Union[str, List, Dict]]
) -> Dict[str, Union[str, Dict]]:
    """Export a workgraph to a rete js editor data."""
    wgdata_short = {
        "name": wgdata["name"],
        "uuid": wgdata.get("uuid", ""),
        "state": wgdata.get("state", ""),
        "nodes": {},
        "links": wgdata.get("links", []),
    }
    #
    for name, task in wgdata["tasks"].items():
        # Add required inputs to nodes
        inputs = []
        for input in task["inputs"]["sockets"].values():
            metadata = input.get("metadata", {}) or {}
            if metadata.get("required", False):
                inputs.append(
                    {"name": input["name"], "identifier": input["identifier"]}
                )

        properties = process_properties(task)
        wgdata_short["nodes"][name] = {
            "label": task["name"],
            "node_type": task["metadata"]["node_type"].upper(),
            "inputs": inputs,
            "properties": properties,
            "outputs": [],
            "position": task["position"],
            "children": task["children"],
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


def validate_task_inout(inout_list: list[str | dict], list_type: str) -> list[dict]:
    """
    Checks if all the list elements provided as `inputs` or `outputs` of to a task are of type `str` or `dict`, and,
    if the former convert them to a list of `dict`s with `name` as the key.

    :param inout_list: The input/output list to be validated.
    :param list_type: "inputs" or "outputs" to indicate what is to be validated for better error message.
    :raises TypeError: If wrong types are provided to the task
    :return: Processed `inputs`/`outputs` list.
    """

    if not all(isinstance(item, (dict, str)) for item in inout_list):
        raise TypeError(
            f"Wrong type provided in the `{list_type}` list to the task, must be either `str` or `dict`."
        )

    processed_inout_list = []

    for item in inout_list:
        if isinstance(item, str):
            processed_inout_list.append({"name": item})
        elif isinstance(item, dict):
            processed_inout_list.append(item)

    processed_inout_list = processed_inout_list

    return processed_inout_list


def filter_keys_namespace_depth(
    dict_: dict[Any, Any], max_depth: int = 0
) -> dict[Any, Any]:
    """
    Filter top-level keys of a dictionary based on the namespace nesting level (number of periods) in the key.

    :param dict dict_: The dictionary to filter.
    :param int max_depth: Maximum depth of namespaces to retain (number of periods).
    :return: The filtered dictionary with only keys satisfying the depth condition.
    :rtype: dict
    """
    result: dict[Any, Any] = {}

    for key, value in dict_.items():
        depth = key.count(".")

        if depth <= max_depth:
            result[key] = value

    return result


def wait_to_link(wgdata: Dict[str, Any]) -> None:
    """Convert wait attribute to link."""
    for name, task in wgdata["tasks"].items():
        for wait_task in task["wait"]:
            if wait_task in wgdata["tasks"]:
                wgdata["links"].append(
                    {
                        "from_node": wait_task,
                        "from_socket": "_wait",
                        "to_node": name,
                        "to_socket": "_wait",
                    }
                )


def shallow_copy_nested_dict(d):
    """Recursively copies only the dictionary structure but keeps value references."""
    if isinstance(d, dict):
        return {key: shallow_copy_nested_dict(value) for key, value in d.items()}
    return d


def make_json_serializable(data):
    """Recursively convert AiiDA objects to JSON-serializable structures."""
    from collections.abc import Mapping, Sequence

    if isinstance(data, orm.Data):
        # Return an int if it's an orm.Int, or a more general dict for other data
        return {
            "__aiida_class__": data.__class__.__name__,
            "uuid": str(data.uuid),
        }
    elif isinstance(data, Mapping):
        return {k: make_json_serializable(v) for k, v in data.items()}
    elif isinstance(data, Sequence) and not isinstance(data, (str, bytes)):
        return [make_json_serializable(item) for item in data]
    else:
        return data
