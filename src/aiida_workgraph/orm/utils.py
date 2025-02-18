from aiida.orm.utils.serialize import _NODE_TAG, _COMPUTER_TAG
from typing import Any
import yaml
from aiida import orm


class AiiDASafeLoader(yaml.SafeLoader):
    """
    A “safe” AiiDA-specific YAML loader.

    Since we are extending SafeLoader, we need to carefully add only those
    constructors we consider safe. Anything we omit will raise an error if
    encountered in the YAML.
    """

    pass


def node_constructor(loader: yaml.SafeLoader, node: yaml.Node) -> orm.Node:
    """
    Load a node from the yaml representation.
    """
    yaml_node = loader.construct_scalar(node)  # type: ignore[arg-type]
    return orm.load_node(uuid=yaml_node)


def computer_constructor(loader: yaml.SafeLoader, computer: yaml.Node) -> orm.Computer:
    """Load a computer from the yaml representation."""
    yaml_node = loader.construct_scalar(computer)  # type: ignore[arg-type]
    return orm.Computer.collection.get(uuid=yaml_node)


AiiDASafeLoader.add_constructor(_NODE_TAG, node_constructor)
AiiDASafeLoader.add_constructor(_COMPUTER_TAG, computer_constructor)


def deserialize_safe(serialized: str) -> Any:

    return yaml.load(serialized, Loader=AiiDASafeLoader)
