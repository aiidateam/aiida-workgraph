from aiida import orm
import importlib.metadata
from typing import Any

builtins_type_mapping = {
    "default": "workgraph.any",
    "namespace": "workgraph.namespace",
    int: "workgraph.int",
    float: "workgraph.float",
    str: "workgraph.string",
    bool: "workgraph.bool",
    list: "workgraph.list",
    dict: "workgraph.dict",
    orm.Int: "workgraph.int",
    orm.Float: "workgraph.float",
    orm.Str: "workgraph.string",
    orm.Bool: "workgraph.bool",
    orm.List: "workgraph.list",
    orm.Dict: "workgraph.dict",
    orm.StructureData: "workgraph.aiida_structuredata",
    Any: "workgraph.any",
}


# Load additional mapping from entry points
def load_custom_type_mapping():
    """Loads custom type mapping from plugins."""
    type_mapping = {}

    entry_points = importlib.metadata.entry_points()

    if hasattr(entry_points, "select"):  # Python 3.10+
        group_entries = entry_points.select(group="aiida_workgraph.type_mapping")
    else:  # Python 3.9 and earlier
        group_entries = entry_points.get("aiida_workgraph.type_mapping", [])

    for entry_point in group_entries:
        try:
            # Load the function or dict and merge with default mapping
            custom_mapping = entry_point.load()
            if isinstance(custom_mapping, dict):
                type_mapping.update(custom_mapping)
        except Exception as e:
            print(f"Failed to load type mapping from {entry_point.name}: {e}")
    return type_mapping


type_mapping = load_custom_type_mapping()
