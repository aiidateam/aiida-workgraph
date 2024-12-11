import json
from aiida.manage import get_config
from pathlib import Path

WORKGRAPH_EXTRA_KEY = "_workgraph"
WORKGRAPH_SHORT_EXTRA_KEY = "_workgraph_short"


builtin_inputs = [
    {"name": "_wait", "link_limit": 1e6, "metadata": {"arg_type": "none"}}
]
builtin_outputs = [{"name": "_wait"}, {"name": "_outputs"}]


def load_config() -> dict:
    """Load the configuration from the config file."""
    config = get_config()
    config_file_path = Path(config.dirpath) / "workgraph.json"
    try:
        with config_file_path.open("r") as f:
            config = json.load(f)
    except FileNotFoundError:
        config = {}
    return config
