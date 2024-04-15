import importlib.metadata
import pathlib

import anywidget
import traitlets
from .utils import wait_to_link

try:
    __version__ = importlib.metadata.version("widget")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"

default_value = {"summary": {}, "nodes": {}, "links": [], "logs": [], "pk": []}


class NodeGraphWidget(anywidget.AnyWidget):
    _esm = pathlib.Path(__file__).parent / "static" / "widget.js"
    _css = pathlib.Path(__file__).parent / "static" / "widget.css"
    value = traitlets.Dict(default_value).tag(sync=True)

    def from_worktree(self, worktree):
        from aiida_worktree.web.backend.app.utils import worktree_to_short_json

        wtdata = worktree.to_dict()
        wait_to_link(wtdata)
        self.value = worktree_to_short_json(wtdata)
