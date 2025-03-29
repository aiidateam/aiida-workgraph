from aiida import orm
import importlib.metadata
from typing import Any

type_mapping = {
    "default": "workgraph.any",
    "namespace": "workgraph.namespace",
    int: "workgraph.int",
    float: "workgraph.float",
    str: "workgraph.string",
    bool: "workgraph.bool",
    orm.Int: "workgraph.aiida_int",
    orm.Float: "workgraph.aiida_float",
    orm.Str: "workgraph.aiida_string",
    orm.Bool: "workgraph.aiida_bool",
    orm.List: "workgraph.aiida_list",
    orm.Dict: "workgraph.aiida_dict",
    orm.StructureData: "workgraph.aiida_structuredata",
    Any: "workgraph.any",
}


# Load additional mapping from entry points
def load_custom_type_mapping():
    """Loads custom type mapping from plugins."""
    for entry_point in importlib.metadata.entry_points().get(
        "workgraph.type_mapping", []
    ):
        try:
            # Load the function or dict and merge with default mapping
            custom_mapping = entry_point.load()
            if isinstance(custom_mapping, dict):
                type_mapping.update(custom_mapping)
        except Exception as e:
            print(f"Failed to load type mapping from {entry_point.name}: {e}")


load_custom_type_mapping()
