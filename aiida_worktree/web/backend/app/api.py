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


@app.get("/worktree/{id}")
async def read_item(id: str):
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
        return {"id": id, "content": content}
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Worktree {id} not found")
