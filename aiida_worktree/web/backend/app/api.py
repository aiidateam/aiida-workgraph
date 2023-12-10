from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from aiida import load_profile
from aiida_worktree.web.backend.app.daemon import router as daemon_router
from aiida_worktree.web.backend.app.worktree import router as worktree_router

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
    return {"message": "Welcome to AiiDA-WorkTree."}


app.include_router(worktree_router)
app.include_router(daemon_router)
