from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from aiida import load_profile
from aiida_worktree.web.backend.app.daemon import router as daemon_router
from aiida_worktree.web.backend.app.worktree import router as worktree_router
from fastapi.staticfiles import StaticFiles
from pathlib import Path
import os

load_profile()

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all methods
    allow_headers=["*"],  # Allows all headers
)


@app.get("/api", tags=["root"])
async def read_root() -> dict:
    return {"message": "Welcome to AiiDA-WorkTree."}


app.include_router(worktree_router)
app.include_router(daemon_router)

# Integrating React build into a FastAPI application and serving the build (HTML, CSS, JavaScript) as static files
# TODO for the moment, one can not refresh the page, because the backend does not know the routes
"""
When you navigate to http://127.0.0.1:8000/settings from http://127.0.0.1:8000/ using client-side
routing (i.e., links within your React app), the React Router handles the route /settings
without reloading the page from the server. This is why it works.
However, when you refresh the page at http://127.0.0.1:8000/settings, the browser makes
a request to the FastAPI server for /settings. Since this route isn't defined in FastAPI
(it's a client-side route), the server returns a 404 Not Found error.
"""
backend_dir = Path(__file__).parent
react_build_directory = backend_dir / "../../frontend/build"
react_build_directory = os.getenv("REACT_BUILD_DIR", react_build_directory)
app.mount("/", StaticFiles(directory=react_build_directory, html=True), name="static")
