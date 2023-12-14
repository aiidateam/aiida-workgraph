from fastapi import APIRouter, HTTPException, Query
from aiida import orm

router = APIRouter()


@router.get("/api/worktree-data")
async def read_worktree_data(search: str = Query(None)):
    from aiida_worktree.cli.query_worktree import WorkTreeQueryBuilder

    try:
        relationships = {}
        builder = WorkTreeQueryBuilder()
        query_set = builder.get_query_set(
            relationships=relationships,
            # Add logic to apply search filter if search query is present
            # This is just an example, you'll need to adjust it based on your actual query builder and data model
            filters={"attributes.process_label": {"like": f"%{search}%"}}
            if search
            else None,
            # order_by, past_days, limit, etc. can also be parameters
        )
        project = ["pk", "uuid", "state", "ctime", "mtime", "process_label"]
        projected = builder.get_projected(query_set, projections=project)
        # pop headers
        projected.pop(0)
        data = []
        for p in projected:
            data.append({project[i]: p[i] for i in range(len(project))})
        return data
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")


@router.get("/api/worktree/{id}/{node_name}")
async def read_worktree_node(id: int, node_name: str):
    from .utils import node_to_short_json
    from aiida.orm.utils.serialize import deserialize_unsafe

    try:
        node = orm.load_node(id)
        wtdata = node.base.extras.get("worktree", None)
        if wtdata is None:
            print("No worktree data found in the node.")
            return

        wtdata = deserialize_unsafe(wtdata)
        content = node_to_short_json(id, wtdata["nodes"][node_name])
        return content
    except KeyError:
        raise HTTPException(
            status_code=404, detail=f"Worktree {id}/{node_name} not found"
        )


@router.get("/api/worktree/{id}")
async def read_worktree_item(id: int):
    from .utils import (
        worktree_to_short_json,
        get_node_summary,
        get_node_inputs,
        get_node_outputs,
    )
    from aiida.cmdline.utils.common import get_workchain_report
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida_worktree.utils import get_parent_worktrees, get_processes_latest

    try:
        node = orm.load_node(id)

        wtdata = node.base.extras.get("worktree", None)
        if wtdata is None:
            print("No worktree data found in the node.")
            return
        wtdata = deserialize_unsafe(wtdata)
        content = worktree_to_short_json(wtdata)
        summary = {
            "table": get_node_summary(node),
            "inputs": get_node_inputs(id),
            "outputs": get_node_outputs(id),
        }
        report = get_workchain_report(node, "REPORT")
        parent_worktrees = get_parent_worktrees(id)
        parent_worktrees.reverse()
        processes_info = get_processes_latest(id)
        content["summary"] = summary
        content["logs"] = report.splitlines()
        content["parent_worktrees"] = parent_worktrees
        content["processes_info"] = processes_info
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")


@router.get("/api/worktree-state/{id}")
async def read_worktree_item_state(id: int):
    from aiida_worktree.utils import get_processes_latest

    try:
        processes_info = get_processes_latest(id)
        return processes_info
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")


# Route for pausing a worktree item
@router.post("/api/worktree/pause/{id}")
async def pause_worktree_node(
    id: int,
):
    from aiida.engine.processes.control import pause_processes

    try:
        # Perform the pause action here
        node = orm.load_node(id)
        pause_processes([node])
        return {"message": f"Send message to pause worktree {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Route for playing a worktree item
@router.post("/api/worktree/play/{id}")
async def play_worktree_node(
    id: int,
):
    from aiida.engine.processes.control import play_processes

    try:
        # Perform the play action here
        node = orm.load_node(id)
        play_processes([node])
        return {"message": f"Send message to play worktree {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Route for deleting a worktree item
@router.delete("/api/worktree/delete/{id}")
async def delete_worktree_node(
    id: int,
):
    from aiida.tools import delete_nodes

    try:
        # Perform the delete action here
        _, was_deleted = delete_nodes([id], dry_run=False)

        if was_deleted:
            return {"message": f"Deleted worktree {id}"}
        else:
            return {"message": f"Failed to delete worktree {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
