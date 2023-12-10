from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from aiida import load_profile, orm

load_profile()

app = FastAPI()

origins = ["http://localhost:3000", "localhost:3000"]


app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/", tags=["root"])
async def read_root() -> dict:
    return {"message": "Welcome to your todo list."}


@app.get("/worktree-data")
async def read_worktree_data():
    from fastapi import HTTPException
    from aiida_worktree.cli.query_worktree import WorkTreeQueryBuilder

    try:
        relationships = {}
        builder = WorkTreeQueryBuilder()
        query_set = builder.get_query_set(
            relationships=relationships,
            # filters=filters,
            # order_by={order_by: order_dir},
            # past_days=past_days,
            # limit=limit,
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


@app.get("/worktree/{id}")
async def read_worktree_item(id: int):
    from fastapi import HTTPException
    from .utils import worktree_to_json

    try:
        node = orm.load_node(id)
        from aiida.orm.utils.serialize import deserialize_unsafe

        wtdata = node.base.extras.get("worktree", None)
        if wtdata is None:
            print("No worktree data found in the node.")
            return
        wtdata = deserialize_unsafe(wtdata)
        content = worktree_to_json(wtdata)
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")
