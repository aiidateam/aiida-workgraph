from typing import Dict, Any


def wait_to_link(wgdata: Dict[str, Any]) -> None:
    """Convert wait attribute to link."""
    for name, node in wgdata["nodes"].items():
        for wait_node in node["wait"]:
            if wait_node in wgdata["nodes"]:
                wgdata["links"].append(
                    {
                        "from_node": wait_node,
                        "from_socket": "_wait",
                        "to_node": name,
                        "to_socket": "_wait",
                    }
                )


def clean_hanging_links(wgdata: Dict[str, Any]) -> None:
    """Clean hanging links in the workgraph."""
    for link in wgdata["links"][:]:  # Iterate over a shallow copy of the list
        if (
            link["from_node"] not in wgdata["nodes"]
            or link["to_node"] not in wgdata["nodes"]
        ):
            wgdata["links"].remove(link)
