from aiida.orm import load_node


def worktree_to_short_json(wtdata):
    """Export a worktree to a rete js editor data."""
    wtdata_short = {
        "name": wtdata["name"],
        "uuid": wtdata["uuid"],
        "state": wtdata["state"],
        "nodes": {},
        "links": wtdata["links"],
    }
    #
    for name, node in wtdata["nodes"].items():
        wtdata_short["nodes"][name] = {
            "label": node["name"],
            "inputs": [],
            "outputs": [],
        }
    for link in wtdata["links"]:
        wtdata_short["nodes"][link["to_node"]]["inputs"].append(
            {
                "name": link["to_socket"],
            }
        )
        wtdata_short["nodes"][link["from_node"]]["outputs"].append(
            {
                "name": link["from_socket"],
            }
        )
    return wtdata_short


def is_function_and_get_source(obj):
    import inspect

    if callable(obj):
        source_lines, _ = inspect.getsourcelines(obj)
        source_code = "".join(source_lines)
        return True, source_code
    else:
        return False, None


def get_node_inputs(pk):
    from aiida.common.links import LinkType
    from aiida.cmdline.utils.common import format_nested_links

    if pk is None:
        return ""

    node = load_node(pk)
    result = ""
    nodes_input = node.base.links.get_incoming(
        link_type=(LinkType.INPUT_CALC, LinkType.INPUT_WORK)
    )
    if nodes_input:
        result += f"{format_nested_links(nodes_input.nested(), headers=['Inputs', 'PK', 'Type'])}"

    return result


def get_node_outputs(pk):
    from aiida.common.links import LinkType
    from aiida.cmdline.utils.common import format_nested_links

    if pk is None:
        return ""

    node = load_node(pk)
    result = ""
    nodes_output = node.base.links.get_outgoing(
        link_type=(LinkType.CREATE, LinkType.RETURN)
    )
    if nodes_output.all():
        result += f"{format_nested_links(nodes_output.nested(), headers=['outputs', 'PK', 'Type'])}"

    return result


def node_to_short_json(worktree_pk, ndata):
    """Export a node to a rete js node."""
    from aiida_worktree.utils import get_executor, get_processes_latest

    executor, _ = get_executor(ndata["executor"])
    is_function, source_code = is_function_and_get_source(executor)
    if is_function:
        executor = source_code
    else:
        executor = str(executor)
    ndata_short = {
        "node_type": ndata["metadata"]["node_type"],
        "metadata": [
            ["name", ndata["name"]],
            ["node_type", ndata["metadata"]["node_type"]],
            ["identifier", ndata["metadata"]["identifier"]],
        ],
        "executor": executor,
    }
    process_info = get_processes_latest(worktree_pk).get(ndata["name"], {})
    ndata_short["process"] = process_info
    if process_info is not None:
        ndata_short["metadata"].append(["pk", process_info.get("pk")])
        ndata_short["metadata"].append(["state", process_info.get("state")])
        ndata_short["metadata"].append(["ctime", process_info.get("ctime")])
        ndata_short["metadata"].append(["mtime", process_info.get("mtime")])
        ndata_short["inputs"] = get_node_inputs(process_info.get("pk"))
        ndata_short["outputs"] = get_node_outputs(process_info.get("pk"))
    else:
        ndata_short["inputs"] = ""
        ndata_short["outputs"] = ""

    return ndata_short


def get_node_summary(node):
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
