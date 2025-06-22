from __future__ import annotations

from typing import Any, Dict, Optional, Union, Callable, List
from aiida.engine.processes import Process
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.engine.runners import Runner
from aiida_workgraph.config import task_types
from aiida.engine import CalcJob, WorkChain
from aiida_pythonjob import PythonJob
from aiida_pythonjob.calculations.pyfunction import PyFunction
from aiida_shell.calculations.shell import ShellJob
import inspect
from yaml.constructor import ConstructorError


def inspect_aiida_component_type(executor: Callable) -> str:
    task_type = None
    if isinstance(executor, type):
        if executor == PythonJob:
            task_type = "PYTHONJOB"
        elif executor == PyFunction:
            task_type = "PYFUNCTION"
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


def get_workgraph_data(
    process: Union[int, orm.Node], safe_load: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Get the workgraph data from the given process node.
    """
    from aiida_workgraph.orm.utils import deserialize_safe
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import load_node

    if safe_load:
        deserializer = deserialize_safe
    else:
        deserializer = deserialize_unsafe

    if isinstance(process, int):
        process = load_node(process)
    wgdata = process.workgraph_data
    task_executors = process.task_executors
    for name, task in wgdata["tasks"].items():
        deserialize_input_values_recursively(task["inputs"], deserializer)
        task["executor"] = task_executors.get(name)
    wgdata["error_handlers"] = process.workgraph_error_handlers
    wgdata["meta_sockets"] = deserializer(wgdata["meta_sockets"])
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
    from aiida_workgraph.orm.workgraph import WorkGraphNode

    tasks = {}
    if pk is None:
        return tasks
    node = aiida.orm.load_node(pk)
    if item_type == "called_process":
        # fetch the process that called by the workgraph
        for link in node.base.links.get_outgoing().all():
            if isinstance(link.node, aiida.orm.ProcessNode):
                tasks[f"{link.link_label}-{link.node.pk}"] = {
                    "pk": link.node.pk,
                    "process_type": link.node.process_type,
                    "state": link.node.process_state.value,
                    "ctime": link.node.ctime,
                    "mtime": link.node.mtime,
                }
    elif item_type == "task":
        if not isinstance(node, WorkGraphNode):
            return tasks
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
    from copy import deepcopy

    wgdata_short = {
        "name": wgdata["name"],
        "uuid": wgdata.get("uuid", ""),
        "state": wgdata.get("state", ""),
        "nodes": {},
        "links": deepcopy(wgdata.get("links", []) + wgdata.get("meta_links", [])),
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
    for name, socket in wgdata["meta_sockets"].items():
        inputs = []
        for input in socket["sockets"].values():
            metadata = input.get("metadata", {}) or {}
            if metadata.get("required", False):
                inputs.append(
                    {"name": input["name"], "identifier": input["identifier"]}
                )
        wgdata_short["nodes"][name] = {
            "label": name,
            "node_type": name,
            "inputs": inputs,
            "properties": {},
            "outputs": [],
            "position": [0, 0],
            "children": [],
        }
    # Add links to nodes
    for link in wgdata_short["links"]:
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
    # hide meta nodes if there is no link to them
    for name, socket in wgdata["meta_sockets"].items():
        node = wgdata_short["nodes"][name]
        if len(node["inputs"]) == 0 and len(node["outputs"]) == 0:
            del wgdata_short["nodes"][name]
    # split ctx node into two nodes, one with only inputs and one with only outputs
    ctx_node = wgdata_short["nodes"].get("graph_ctx")
    if ctx_node:
        ctx_inputs = []
        ctx_outputs = []
        for input in ctx_node["inputs"]:
            if input.get("identifier") == "_wait":
                continue
            ctx_inputs.append(input)
        for output in ctx_node["outputs"]:
            if output.get("identifier") == "_wait":
                continue
            ctx_outputs.append(output)
        wgdata_short["nodes"]["ctx_inputs"] = {
            "label": "ctx_inputs",
            "node_type": "graph_ctx",
            "inputs": ctx_inputs,
            "properties": {},
            "outputs": [],
            "position": [0, 0],
            "children": [],
        }
        wgdata_short["nodes"]["ctx_outputs"] = {
            "label": "ctx_outputs",
            "node_type": "graph_ctx",
            "inputs": [],
            "properties": {},
            "outputs": ctx_outputs,
            "position": [0, 0],
            "children": [],
        }
        del wgdata_short["nodes"]["graph_ctx"]
        # update the links
        for link in wgdata_short["links"]:
            if link["from_node"] == "graph_ctx":
                link["from_node"] = "ctx_inputs"
            if link["to_node"] == "graph_ctx":
                link["to_node"] = "ctx_outputs"
    return wgdata_short


def remove_output_values(outputs: Dict[str, Any]) -> None:
    """
    Remove output values from the outputs dictionary.

    This function iterates through the outputs and removes the 'value' key from each output.
    It is useful for cleaning up outputs before serialization or storage.

    :param outputs: A dictionary of outputs to be cleaned.
    :return: None
    """
    if "property" in outputs:
        outputs["property"].pop("value", None)
    if "sockets" in outputs:
        for socket in outputs["sockets"].values():
            remove_output_values(socket)


def serialize_input_values_recursively(
    inputs: Dict[str, Any], serializer: callable = None
) -> None:
    """
    Serialize input values to a format suitable for storage or transmission.

    This function converts all input values to their raw representations, ensuring
    that complex objects are converted to simple types (e.g., strings, integers).

    :param inputs: A dictionary of inputs to be serialized.
    :return: A dictionary with serialized input values.
    """
    if serializer is None:
        from aiida.orm.utils.serialize import serialize

        serializer = serialize
    if "property" in inputs:
        inputs["property"]["value"] = serializer(inputs["property"]["value"])
    if "sockets" in inputs:
        for socket in inputs["sockets"].values():
            serialize_input_values_recursively(socket)


def deserialize_input_values_recursively(
    inputs: Dict[str, Any], deserializer: callable = None
) -> None:
    """
    Deserialize input values from a format suitable for storage or transmission.

    This function converts all input values back to their original representations,
    ensuring that complex objects are restored from their serialized forms.

    :param inputs: A dictionary of inputs to be deserialized.
    :return: A dictionary with deserialized input values.
    """
    if deserializer is None:
        from aiida_workgraph.orm.utils import deserialize_safe

        deserializer = deserialize_safe

    if "property" in inputs:
        try:
            inputs["property"]["value"] = deserializer(inputs["property"]["value"])
        except ConstructorError:
            name = inputs["name"]
            print(
                f"Info: could not deserialize input '{name}'. The raw input value is left in place."
                "The workgraph is still loaded and you can inspect tasks, inputs and outputs. "
                "If you trust the data, reload with safe_load=False."
            )
    if "sockets" in inputs:
        for socket in inputs["sockets"].values():
            deserialize_input_values_recursively(socket, deserializer)


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


def query_existing_processes(pks: list[int]) -> list[int]:
    """Query all existing processes from the database."""
    from aiida import orm

    qb = orm.QueryBuilder()
    qb.append(
        orm.ProcessNode,
        filters={"id": {"in": pks}},
        project=["id"],
    )
    results = qb.all()
    existing_pks = [res[0] for res in results]
    return existing_pks


def query_terminated_processes(pks: list[int]) -> list[int]:
    """Query all terminated processes from the database."""
    from aiida import orm

    qb = orm.QueryBuilder()
    qb.append(
        orm.ProcessNode,
        filters={
            "id": {"in": pks},
            "attributes.process_state": {"in": ["killed", "finished", "excepted"]},
        },
        project=["id"],
    )
    results = qb.all()
    terminated_pks = [res[0] for res in results]
    return terminated_pks
