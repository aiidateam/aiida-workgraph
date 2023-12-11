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


@router.get("/api/worktree/{id}")
async def read_worktree_item(id: int):
    from .utils import worktree_to_short_json

    try:
        node = orm.load_node(id)
        from aiida.orm.utils.serialize import deserialize_unsafe

        wtdata = node.base.extras.get("worktree", None)
        if wtdata is None:
            print("No worktree data found in the node.")
            return
        wtdata = deserialize_unsafe(wtdata)
        content = worktree_to_short_json(wtdata)
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")
