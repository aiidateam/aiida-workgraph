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
    states = traitlets.Dict({}).tag(sync=True)
    positions = traitlets.Dict({}).tag(sync=True)

    def __init__(self, parent=None, **kwargs):
        self.parent = parent
        super().__init__(**kwargs)

    def from_workgraph(self, workgraph):
        from aiida_workgraph.web.backend.app.utils import workgraph_to_short_json

        wgdata = workgraph.to_dict()
        wait_to_link(wgdata)
        wgdata = workgraph_to_short_json(wgdata)
        self.value = wgdata

    def from_node(self, node):
        ndata = node.to_dict()
        ndata.pop("properties", None)
        ndata["label"] = ndata["metadata"]["identifier"]
        wgdata = {"nodes": {node.name: ndata}, "links": []}
        self.value = wgdata

    @traitlets.observe("positions")
    def _observe_positions(self, change):
        if not self.parent:
            return
        if change["new"]:
            for name, pos in change["new"].items():
                self.parent.nodes[name].position = pos
