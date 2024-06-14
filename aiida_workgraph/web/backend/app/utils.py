from typing import Dict, Optional, Union, Tuple, List, Any
from aiida.orm import load_node, Node
from datetime import datetime
from dateutil import relativedelta
from dateutil.tz import tzlocal


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


def is_function_and_get_source(obj: Any) -> Tuple[bool, Optional[str]]:
    import inspect

    if callable(obj):
        source_lines, _ = inspect.getsourcelines(obj)
        source_code = "".join(source_lines)
        return True, source_code
    else:
        return False, None


def get_node_recursive(links: Dict) -> Dict[str, Union[List[int], str]]:
    """Recursively get a dictionary of nodess."""
    from collections.abc import Mapping

    data = {}
    for label, value in links.items():
        if isinstance(value, Mapping):
            data.update({label: get_node_recursive(value)})
        else:
            data[label] = [value.pk, value.__class__.__name__]
    return data


def get_node_inputs(pk: Optional[int]) -> Union[str, Dict[str, Union[List[int], str]]]:
    from aiida.common.links import LinkType

    if pk is None:
        return ""

    node = load_node(pk)
    nodes_input = node.base.links.get_incoming(
        link_type=(LinkType.INPUT_CALC, LinkType.INPUT_WORK)
    )
    if nodes_input:
        result = get_node_recursive(nodes_input.nested())
    else:
        result = {}

    return result


def get_node_outputs(pk: Optional[int]) -> Union[str, Dict[str, Union[List[int], str]]]:
    from aiida.common.links import LinkType

    if pk is None:
        return ""

    node = load_node(pk)
    result = ""
    nodes_output = node.base.links.get_outgoing(
        link_type=(LinkType.CREATE, LinkType.RETURN)
    )
    if nodes_output.all():
        result = get_node_recursive(nodes_output.nested())
    else:
        result = {}

    return result


def node_to_short_json(workgraph_pk: int, tdata: Dict[str, Any]) -> Dict[str, Any]:
    """Export a node to a rete js node."""
    from aiida_workgraph.utils import get_executor, get_processes_latest

    executor, _ = get_executor(tdata["executor"])
    is_function, source_code = is_function_and_get_source(executor)
    if is_function:
        executor = source_code
    else:
        executor = str(executor)
    tdata_short = {
        "node_type": tdata["metadata"]["node_type"],
        "metadata": [
            ["name", tdata["name"]],
            ["node_type", tdata["metadata"]["node_type"]],
            ["identifier", tdata["metadata"]["identifier"]],
        ],
        "executor": executor,
    }
    process_info = get_processes_latest(workgraph_pk).get(tdata["name"], {})
    tdata_short["process"] = process_info
    if process_info is not None:
        tdata_short["metadata"].append(["pk", process_info.get("pk")])
        tdata_short["metadata"].append(["state", process_info.get("state")])
        tdata_short["metadata"].append(["ctime", process_info.get("ctime")])
        tdata_short["metadata"].append(["mtime", process_info.get("mtime")])
        tdata_short["inputs"] = get_node_inputs(process_info.get("pk"))
        tdata_short["outputs"] = get_node_outputs(process_info.get("pk"))
    else:
        tdata_short["inputs"] = ""
        tdata_short["outputs"] = ""

    return tdata_short


def get_node_summary(node: Node) -> List[List[str]]:
    """ """
    from plumpy import ProcessState
    from aiida.orm import ProcessNode

    table = []

    if isinstance(node, ProcessNode):
        table.append(["type", node.process_label])

        try:
            process_state = ProcessState(node.process_state)
        except (AttributeError, ValueError):
            pass
        else:
            process_state_string = process_state.value.capitalize()

            if process_state == ProcessState.FINISHED and node.exit_message:
                table.append(
                    [
                        "state",
                        f"{process_state_string} [{node.exit_status}] {node.exit_message}",
                    ]
                )
            elif process_state == ProcessState.FINISHED:
                table.append(["state", f"{process_state_string} [{node.exit_status}]"])
            elif process_state == ProcessState.EXCEPTED:
                table.append(["state", f"{process_state_string} <{node.exception}>"])
            else:
                table.append(["state", process_state_string])

    else:
        table.append(["type", node.__class__.__name__])
    table.append(["pk", str(node.pk)])
    table.append(["uuid", str(node.uuid)])
    table.append(["label", node.label])
    table.append(["description", node.description])
    table.append(["ctime", node.ctime])
    table.append(["mtime", node.mtime])

    try:
        computer = node.computer
    except AttributeError:
        pass
    else:
        if computer is not None:
            table.append(["computer", f"[{node.computer.pk}] {node.computer.label}"])

    return table


def time_ago(past_time: datetime) -> str:
    # Get the current time
    now = datetime.now(tzlocal())

    # Calculate the time difference
    delta = relativedelta.relativedelta(now, past_time)

    # Format the time difference
    if delta.years > 0:
        return f"{delta.years}Y ago"
    elif delta.months > 0:
        return f"{delta.months}M ago"
    elif delta.days > 0:
        return f"{delta.days}D ago"
    elif delta.hours > 0:
        return f"{delta.hours}h ago"
    elif delta.minutes > 0:
        return f"{delta.minutes}min ago"
    else:
        return "Just now"
