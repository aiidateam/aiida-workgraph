from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from aiida import load_profile
from aiida_workgraph.web.backend.app.daemon import router as daemon_router
from aiida_workgraph.web.backend.app.workgraph import router as workgraph_router
from aiida_workgraph.web.backend.app.datanode import router as datanode_router
from fastapi.staticfiles import StaticFiles
from pathlib import Path
import os

from fastapi.responses import FileResponse
from fastapi.exception_handlers import http_exception_handler
from starlette.exceptions import HTTPException as StarletteHTTPException


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
    return {"message": "Welcome to AiiDA-WorkGraph."}


app.include_router(workgraph_router)
app.include_router(datanode_router)
app.include_router(daemon_router)

# Integrating React build into a FastAPI application and serving the build (HTML, CSS, JavaScript) as static files
"""
When navigate to http://127.0.0.1:8000/settings from http://127.0.0.1:8000/ using client-side
routing (i.e., links within your React app), the React Router handles the route /settings
without reloading the page from the server. This is why it works.
However, when you refresh the page at http://127.0.0.1:8000/settings, the browser makes
a request to the FastAPI server for /settings. Since this route isn't defined in FastAPI
(it's a client-side route), the server returns a 404 Not Found error.
so we use the index.html serve all routes except API specific ones, then load all static assets.
"""
backend_dir = Path(__file__).parent
build_dir = backend_dir / "../../frontend/build/"
build_dir = os.getenv("REACT_BUILD_DIR", build_dir)


@app.exception_handler(StarletteHTTPException)
async def _spa_server(req: Request, exc: StarletteHTTPException):
    if exc.status_code == 404:
        return FileResponse(f"{build_dir}/index.html", media_type="text/html")
    else:
        return await http_exception_handler(req, exc)


if os.path.isdir(build_dir):
    app.mount(
        "/static/",
        StaticFiles(directory=build_dir / "static"),
        name="React app static files",
    )
