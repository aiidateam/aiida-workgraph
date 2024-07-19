from .general_data import GeneralData
from aiida import orm, common
from importlib.metadata import entry_points
from typing import Any
from aiida_workgraph.config import load_config
import sys


def get_serializer_from_entry_points() -> dict:
    """Retrieve the serializer from the entry points."""
    # import time

    # ts = time.time()
    configs = load_config()
    serializers = configs.get("serializers", {})
    excludes = serializers.get("excludes", [])
    # Retrieve the entry points for 'aiida.data' and store them in a dictionary
    eps = entry_points()
    if sys.version_info >= (3, 10):
        group = eps.select(group="aiida.data")
    else:
        group = eps.get("aiida.data", [])
    eps = {}
    for ep in group:
        # split the entry point name by first ".", and check the last part
        key = ep.name.split(".", 1)[-1]
        # skip key without "." because it is not a module name for a data type
        if "." not in key or key in excludes:
            continue
        eps.setdefault(key, [])
        eps[key].append(ep)

    # print("Time to load entry points: ", time.time() - ts)
    # check if there are duplicates
    for key, value in eps.items():
        if len(value) > 1:
            if key in serializers:
                [ep for ep in value if ep.name == serializers[key]]
                eps[key] = [ep for ep in value if ep.name == serializers[key]]
                if not eps[key]:
                    raise ValueError(
                        f"Entry point {serializers[key]} not found for {key}"
                    )
            else:
                msg = f"Duplicate entry points for {key}: {[ep.name for ep in value]}"
                raise ValueError(msg)
    return eps


eps = get_serializer_from_entry_points()


def serialize_to_aiida_nodes(inputs: dict = None) -> dict:
    """Serialize the inputs to a dictionary of AiiDA data nodes.

    Args:
        inputs (dict): The inputs to be serialized.

    Returns:
        dict: The serialized inputs.
    """
    new_inputs = {}
    # save all kwargs to inputs port
    for key, data in inputs.items():
        new_inputs[key] = general_serializer(data)
    return new_inputs


def clean_dict_key(data):
    """Replace "." with "__dot__" in the keys of a dictionary."""
    if isinstance(data, dict):
        return {k.replace(".", "__dot__"): clean_dict_key(v) for k, v in data.items()}
    return data


def general_serializer(data: Any, check_value=True) -> orm.Node:
    """Serialize the data to an AiiDA data node."""
    if isinstance(data, orm.Data):
        if check_value and not hasattr(data, "value"):
            raise ValueError("Only AiiDA data Node with a value attribute is allowed.")
        return data
    elif isinstance(data, common.extendeddicts.AttributeDict):
        # if the data is an AttributeDict, use it directly
        return data
    # if is string with syntax {{}}, this is a port will read data from ctx
    elif isinstance(data, str) and data.startswith("{{") and data.endswith("}}"):
        return data
    # if data is a class instance, get its __module__ and class name as a string
    # for example, an Atoms will have ase.atoms.Atoms
    else:
        data = clean_dict_key(data)
        # try to get the serializer from the entry points
        data_type = type(data)
        ep_key = f"{data_type.__module__}.{data_type.__name__}"
        # search for the key in the entry points
        if ep_key in eps:
            try:
                new_node = eps[ep_key][0].load()(data)
            except Exception as e:
                raise ValueError(f"Error in serializing {ep_key}: {e}")
            finally:
                # try to save the node to da
                try:
                    new_node.store()
                    return new_node
                except Exception:
                    # try to serialize the value as a GeneralData
                    try:
                        new_node = GeneralData(data)
                        new_node.store()
                        return new_node
                    except Exception as e:
                        raise ValueError(f"Error in serializing {ep_key}: {e}")
        else:
            # try to serialize the data as a GeneralData
            try:
                new_node = GeneralData(data)
                new_node.store()
                return new_node
            except Exception as e:
                raise ValueError(f"Error in serializing {ep_key}: {e}")
