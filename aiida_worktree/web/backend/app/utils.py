def worktree_to_short_json(wtdata):
    """Export a worktree to a rete js editor data."""
    wtdata_short = {
        "name": wtdata["name"],
        "uuid": wtdata["uuid"],
        "state": wtdata["state"],
        "nodes": {},
        "links": wtdata["links"],
    }
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
