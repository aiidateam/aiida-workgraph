from typing import Optional, Dict, Any
import importlib.metadata
import pathlib

import anywidget
import traitlets
from .utils import wait_to_link
from aiida_workgraph.utils import filter_keys_namespace_depth

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

    def __init__(self, parent: Optional[Any] = None, **kwargs: Any) -> None:
        self.parent = parent
        super().__init__(**kwargs)

    def from_workgraph(self, workgraph: Any) -> None:
        from aiida_workgraph.utils import workgraph_to_short_json

        wgdata = workgraph.to_dict()
        wait_to_link(wgdata)
        wgdata = workgraph_to_short_json(wgdata)
        self.value = wgdata

    def from_node(self, node: Any, show_socket_depth: int = 0) -> None:

        tdata = node.to_dict()

        # Remove certain elements of the dict-representation of the Node that we don't want to show
        for key in ("properties", "executor", "node_class", "process"):
            tdata.pop(key, None)
        for input in tdata["inputs"].values():
            input.pop("property")

        tdata["label"] = tdata["identifier"]

        filtered_inputs = filter_keys_namespace_depth(
            dict_=tdata["inputs"], max_depth=show_socket_depth
        )
        tdata["inputs"] = list(filtered_inputs.values())
        tdata["outputs"] = list(tdata["outputs"].values())
        wgdata = {"name": node.name, "nodes": {node.name: tdata}, "links": []}
        self.value = wgdata

    def to_html(self, output: str = None, width: str = "100%", height: str = "600px"):
        """Write a standalone html file to visualize the workgraph."""
        from IPython.display import IFrame
        from .html_template import html_template
        import json

        if output is None:
            # create "html" folder if it does not exist
            pathlib.Path("html").mkdir(exist_ok=True)
            output = f"html/{self.value['name']}.html"
        # Replace the placeholder with the actual workgraphData
        html_content = html_template.replace(
            "__WORKGRAPH_DATA__", json.dumps(self.value)
        )
        with open(output, "w") as f:
            f.write(html_content)
        return IFrame(output, width=width, height=height)

    @traitlets.observe("positions")
    def _observe_positions(self, change: Dict[str, Any]) -> None:
        if not self.parent:
            return
        if change["new"]:
            for name, pos in change["new"].items():
                self.parent.nodes[name].position = pos
