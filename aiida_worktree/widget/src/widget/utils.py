def wait_to_link(wtdata):
    """Convert wait attribute to link."""
    for name, node in wtdata["nodes"].items():
        for wait_node in node["wait"]:
            if wait_node in wtdata["nodes"]:
                wtdata["links"].append(
                    {
                        "from_node": wait_node,
                        "from_socket": "_wait",
                        "to_node": name,
                        "to_socket": "_wait",
                    }
                )


def clean_hanging_links(wtdata):
    """Clean hanging links in the nodetree."""
    for link in wtdata["links"][:]:  # Iterate over a shallow copy of the list
        if (
            link["from_node"] not in wtdata["nodes"]
            or link["to_node"] not in wtdata["nodes"]
        ):
            wtdata["links"].remove(link)
