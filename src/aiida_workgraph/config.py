import json
from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER

WORKGRAPH_EXTRA_KEY = "_workgraph"
WORKGRAPH_SHORT_EXTRA_KEY = "_workgraph_short"


builtin_inputs = [{"name": "_wait", "link_limit": 1e6, "arg_type": "none"}]
builtin_outputs = [{"name": "_wait"}, {"name": "_outputs"}]


def load_config() -> dict:
    """Load the configuration from the config file."""
    config_file_path = AIIDA_CONFIG_FOLDER / "workgraph.json"
    try:
        with config_file_path.open("r") as f:
            config = json.load(f)
    except FileNotFoundError:
        config = {}
    return config
