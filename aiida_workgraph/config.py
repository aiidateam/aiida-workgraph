import json
from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER


def load_config() -> dict:
    """Load the configuration from the config file."""
    config_file_path = AIIDA_CONFIG_FOLDER / "workgraph.json"
    try:
        with config_file_path.open("r") as f:
            config = json.load(f)
    except FileNotFoundError:
        config = {}
    return config
