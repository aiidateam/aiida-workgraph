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
default_style = {"width": "90%", "height": "600px"}
default_settings = {"minimap": True}


class NodeGraphWidget(anywidget.AnyWidget):
    _esm = pathlib.Path(__file__).parent / "static" / "widget.js"
    _css = pathlib.Path(__file__).parent / "static" / "widget.css"
    value = traitlets.Dict(default_value).tag(sync=True)
    settings = traitlets.Dict(default_settings).tag(sync=True)
    style = traitlets.Dict(default_style).tag(sync=True)

    def from_worktree(self, worktree):
        from aiida_worktree.web.backend.app.utils import worktree_to_short_json

        wtdata = worktree.to_dict()
        wait_to_link(wtdata)
        wtdata = worktree_to_short_json(wtdata)
        self.value = wtdata

    def from_node(self, node):
        ndata = node.to_dict()
        ndata["label"] = ndata["metadata"]["identifier"]
        wtdata = {"nodes": {node.name: ndata}, "links": []}
        self.value = wtdata
