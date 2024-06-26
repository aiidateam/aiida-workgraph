from fastapi import APIRouter, HTTPException, Query
from aiida import orm
from typing import List

router = APIRouter()


@router.get("/api/workgraph-data")
async def read_workgraph_data(search: str = Query(None)):
    from aiida_workgraph.cli.query_workgraph import WorkGraphQueryBuilder

    try:
        relationships = {}
        builder = WorkGraphQueryBuilder()
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
        raise HTTPException(status_code=404, detail=f"Workgraph {id} not found")


@router.get("/api/workgraph/{id}/{node_name}")
async def read_workgraph_task(id: int, node_name: str):
    from .utils import node_to_short_json
    from aiida.orm.utils.serialize import deserialize_unsafe

    try:
        node = orm.load_node(id)
        wgdata = node.base.extras.get("_workgraph", None)
        if wgdata is None:
            print("No workgraph data found in the node.")
            return

        wgdata = deserialize_unsafe(wgdata)
        content = node_to_short_json(id, wgdata["tasks"][node_name])
        return content
    except KeyError:
        raise HTTPException(
            status_code=404, detail=f"Workgraph {id}/{node_name} not found"
        )


@router.get("/api/workgraph/{id}")
async def read_workgraph(id: int):
    from .utils import (
        workgraph_to_short_json,
        get_node_summary,
        get_node_inputs,
        get_node_outputs,
    )
    from aiida.orm.utils.serialize import deserialize_unsafe
    from aiida_workgraph.utils import get_parent_workgraphs, get_processes_latest

    try:
        node = orm.load_node(id)

        wgdata = node.base.extras.get("_workgraph", None)
        if wgdata is None:
            print("No workgraph data found in the node.")
            return
        wgdata = deserialize_unsafe(wgdata)
        content = workgraph_to_short_json(wgdata)
        summary = {
            "table": get_node_summary(node),
            "inputs": get_node_inputs(id),
            "outputs": get_node_outputs(id),
        }

        parent_workgraphs = get_parent_workgraphs(id)
        parent_workgraphs.reverse()
        processes_info = get_processes_latest(id)
        content["summary"] = summary
        content["parent_workgraphs"] = parent_workgraphs
        content["processes_info"] = processes_info
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Workgraph {id} not found")


@router.get("/api/workgraph-state/{id}")
async def read_workgraph_tasks_state(id: int):
    from aiida_workgraph.utils import get_processes_latest

    try:
        processes_info = get_processes_latest(id)
        return processes_info
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Workgraph {id} not found")


@router.get("/api/workgraph-logs/{id}")
async def read_workgraph_logs(id: int):
    from aiida.cmdline.utils.common import get_workchain_report

    try:
        node = orm.load_node(id)
        report = get_workchain_report(node, "REPORT")
        logs = report.splitlines()
        return logs
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Workgraph {id} not found")


# Route for pausing a workgraph item
@router.post("/api/workgraph/pause/{id}")
async def pause_workgraph(
    id: int,
):
    from aiida.engine.processes.control import pause_processes

    try:
        # Perform the pause action here
        node = orm.load_node(id)
        pause_processes([node])
        return {"message": f"Send message to pause workgraph {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Route for playing a workgraph item
@router.post("/api/workgraph/play/{id}")
async def play_workgraph(
    id: int,
):
    from aiida.engine.processes.control import play_processes

    try:
        # Perform the play action here
        node = orm.load_node(id)
        play_processes([node])
        return {"message": f"Send message to play workgraph {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Route for deleting a workgraph item
@router.delete("/api/workgraph/delete/{id}")
async def delete_workgraph(
    id: int,
):
    from aiida.tools import delete_nodes

    try:
        # Perform the delete action here
        _, was_deleted = delete_nodes([id], dry_run=False)

        if was_deleted:
            return {"message": f"Deleted workgraph {id}"}
        else:
            return {"message": f"Failed to delete workgraph {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# General function to manage task actions
async def manage_task_action(action: str, id: int, tasks: List[str]):
    from aiida_workgraph.utils.control import pause_tasks, play_tasks, kill_tasks

    print(f"Performing {action} action on tasks {tasks} in workgraph {id}")
    try:

        if action == "pause":
            print(f"Pausing tasks {tasks}")
            _, msg = pause_tasks(id, tasks=tasks)
        elif action == "play":
            print(f"Playing tasks {tasks}")
            _, msg = play_tasks(id, tasks)
        elif action == "kill":
            print(f"Killing tasks {tasks}")
            _, msg = kill_tasks(id, tasks)
        else:
            raise HTTPException(status_code=400, detail="Unsupported action")

        return {"message": msg}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Endpoint for pausing tasks in a workgraph
@router.post("/api/workgraph/tasks/pause/{id}")
async def pause_workgraph_tasks(id: int, tasks: List[str] = None):
    return await manage_task_action("pause", id, tasks)


# Endpoint for playing tasks in a workgraph
@router.post("/api/workgraph/tasks/play/{id}")
async def play_workgraph_tasks(id: int, tasks: List[str] = None):
    return await manage_task_action("play", id, tasks)


# Endpoint for killing tasks in a workgraph
@router.post("/api/workgraph/tasks/kill/{id}")
async def kill_workgraph_tasks(id: int, tasks: List[str] = None):
    return await manage_task_action("kill", id, tasks)
