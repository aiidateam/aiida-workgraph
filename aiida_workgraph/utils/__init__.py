from __future__ import annotations

from typing import Any, Dict, Optional, Union, Callable, List
from aiida.engine.processes import Process
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.engine.runners import Runner


def build_callable(obj: Callable) -> Dict[str, Any]:
    """
    Build the executor data from the callable. This will either serialize the callable
    using cloudpickle if it's a local or lambda function, or store its module and name
    if it's a globally defined callable (function or class).

    Args:
        obj (Callable): The callable to be serialized or referenced.

    Returns:
        Dict[str, Any]: A dictionary containing the serialized callable or a reference
                        to its module and name.
    """
    import types
    from aiida_workgraph.orm.function_data import PickledFunction

    # Check if the callable is a function or class
    if isinstance(obj, (types.FunctionType, type)):
        # Check if callable is nested (contains dots in __qualname__ after the first segment)
        if obj.__module__ == "__main__" or "." in obj.__qualname__.split(".", 1)[-1]:
            # Local or nested callable, so pickle the callable
            executor = PickledFunction.build_callable(obj)
        else:
            # Global callable (function/class), store its module and name for reference
            executor = {
                "module": obj.__module__,
                "name": obj.__name__,
                "is_pickle": False,
            }
    elif isinstance(obj, PickledFunction) or isinstance(obj, dict):
        executor = obj
    else:
        raise TypeError("Provided object is not a callable function or class.")
    return executor


def get_sorted_names(data: dict) -> list[str]:
    """Get the sorted names from a dictionary."""
    sorted_names = [
        name
        for name, _ in sorted(
            ((name, item["list_index"]) for name, item in data.items()),
            key=lambda x: x[1],
        )
    ]
    return sorted_names


def store_nodes_recursely(data: Any) -> None:
    """Recurse through a data structure and store any unstored nodes that are found along the way
    :param data: a data structure potentially containing unstored nodes
    """
    from aiida.orm import Node
    import collections.abc

    if isinstance(data, Node) and not data.is_stored:
        data.store()
    elif isinstance(data, collections.abc.Mapping):
        for _, value in data.items():
            store_nodes_recursely(value)
    elif isinstance(data, collections.abc.Sequence) and not isinstance(data, str):
        for value in data:
            store_nodes_recursely(value)


def get_executor(data: Dict[str, Any]) -> Union[Process, Any]:
    """Import executor from path and return the executor and type."""
    import importlib
    from aiida.plugins import CalculationFactory, WorkflowFactory, DataFactory

    data = data or {}
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
        elif not data.get("name", None) and not data.get("module", None):
            executor = None
        else:
            module = importlib.import_module("{}".format(data.get("module", "")))
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


def organize_nested_inputs(wgdata: Dict[str, Any]) -> None:
    """Merge sub properties to the root properties.
    The sub properties will be se
    For example:
        task["inputs"]["base"]["property"]["value"] = None
        task["inputs"]["base.pw.parameters"]["property"]["value"] = 2
        task["inputs"]["base.pw.code"]["property"]["value"] = 1
        task["inputs"]["metadata"]["property"]["value"] = {"options": {"resources": {"num_cpus": 1}}
        task["inputs"]["metadata.options"]["property"]["value"] = {"resources": {"num_machine": 1}}
    After organizing:
        task["inputs"]["base"]["property"]["value"] = {"base": {"pw": {"parameters": 2,
                                                                     "code": 1},
                                                       "metadata": {"options":
                                                                        {"resources": {"num_cpus": 1,
                                                                                       "num_machine": 1}}}},
                                                       }
        task["inputs"]["base.pw.parameters"]["property"]["value"] = None
        task["inputs"]["base.pw.code"]["property"]["value"] = None
        task["inputs"]["metadata"]["property"]["value"] = None
        task["inputs"]["metadata.options"]["property"]["value"] = None
    """
    for _, task in wgdata["tasks"].items():
        for key, prop in task["properties"].items():
            if "." in key and prop["value"] not in [None, {}]:
                root, key = key.split(".", 1)
                root_prop = task["properties"][root]
                update_nested_dict(root_prop["value"], key, prop["value"])
                prop["value"] = None
        for key, input in task["inputs"].items():
            if input["property"] is None:
                continue
            prop = input["property"]
            if "." in key and prop["value"] not in [None, {}]:
                root, key = key.split(".", 1)
                root_prop = task["inputs"][root]["property"]
                # update the root property
                root_prop["value"] = update_nested_dict(
                    root_prop["value"], key, prop["value"]
                )
                prop["value"] = None


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
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import load_node

    if isinstance(process, int):
        process = load_node(process)
    wgdata = process.base.extras.get("_workgraph", None)
    if wgdata is None:
        return
    for name, task in wgdata["tasks"].items():
        wgdata["tasks"][name] = deserialize_unsafe(task)
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
    if len(links) > 0:
        parent_workgraphs.extend(get_parent_workgraphs(links[0].node.pk))
    return parent_workgraphs


def get_processes_latest(
    pk: int, node_name: str = None, item_type: str = "task"
) -> Dict[str, Dict[str, Union[int, str]]]:
    """Get the latest info of all tasks from the process."""
    import aiida
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida.orm import QueryBuilder
    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    tasks = {}
    if item_type == "called_process":
        # fetch the process that called by the workgraph
        node = aiida.orm.load_node(pk)
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


def serialize_properties(wgdata):
    """Serialize the properties.
    Because we use yaml (aiida's serialize) to serialize the data and
    save it to the node.base.extras. yaml can not handle the function
    defined in a scope, e.g., local function in another function.
    So, if a function is used as input, we needt to serialize the function.
    """
    from aiida_workgraph.orm.function_data import PickledLocalFunction
    from aiida_workgraph.tasks.pythonjob import PythonJob
    import inspect

    for _, task in wgdata["tasks"].items():
        if task["metadata"]["node_type"].upper() == "PYTHONJOB":
            PythonJob.serialize_pythonjob_data(task)
        for _, input in task["inputs"].items():
            if input["property"] is None:
                continue
            prop = input["property"]
            if inspect.isfunction(prop["value"]):
                prop["value"] = PickledLocalFunction(prop["value"]).store()


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
    for name, input in task["inputs"].items():
        if input["property"] is not None:
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
        "uuid": wgdata["uuid"],
        "state": wgdata["state"],
        "nodes": {},
        "links": wgdata["links"],
    }
    #
    for name, task in wgdata["tasks"].items():
        # Add required inputs to nodes
        inputs = [
            {"name": name, "identifier": input["identifier"]}
            for name, input in task["inputs"].items()
            if name in task["args"]
            or (task["identifier"].upper() == "SHELLJOB" and name.startswith("nodes."))
        ]

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
