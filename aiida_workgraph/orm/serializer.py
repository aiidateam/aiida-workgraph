from .general_data import GeneralData
from aiida import orm, common
from importlib.metadata import entry_points
from typing import Any


# Retrieve the entry points for 'aiida.data' and store them in a dictionary
eps = {ep.name: ep for ep in entry_points().get("aiida.data", [])}


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
                new_node = eps[ep_key].load()(data)
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
