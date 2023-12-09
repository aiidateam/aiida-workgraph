def worktree_to_json(wtdata):
    """Export a worktree to a JSON serializable dictionary."""
    wtdata_short = {
        "name": wtdata["name"],
        "uuid": wtdata["uuid"],
        "state": wtdata["state"],
        "nodes": {},
        "links": wtdata["links"],
    }
    for name, node in wtdata["nodes"].items():
        wtdata_short["nodes"][name] = {
            "name": node["name"],
            "identifier": node["metadata"]["identifier"],
            "node_type": node["metadata"]["node_type"],
            "uuid": node["uuid"],
            "state": node["state"],
            "action": node["action"],
            "inputs": node["inputs"],
        }
    return wtdata_short
