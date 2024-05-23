from .general_data import GeneralData
from aiida import orm
from importlib.metadata import entry_points


# Retrieve the entry points for 'aiida.data' and store them in a dictionary
eps = {ep.name: ep for ep in entry_points().get("aiida.data", [])}


def general_serializer(inputs):
    """Serialize the inputs to a dictionary of AiiDA data nodes.

    Args:
        inputs (dict): The inputs to be serialized.

    Returns:
        dict: The serialized inputs.
    """
    new_inputs = {}
    # save all kwargs to inputs port
    for key, value in inputs.items():
        if isinstance(value, orm.Data):
            if not hasattr(value, "value"):
                raise ValueError(
                    "Only AiiDA data Node with a value attribute is allowed."
                )
            new_inputs[key] = value
        # if value is a class instance, get its __module__ and class name as a string
        # for example, an Atoms will have ase.atoms.Atoms
        else:
            # try to get the serializer from the entry points
            value_type = type(value)
            ep_key = f"{value_type.__module__}.{value_type.__name__}"
            # search for the key in the entry points
            if ep_key in eps:
                try:
                    new_inputs[key] = eps[ep_key].load()(value)
                except Exception as e:
                    raise ValueError(f"Error in serializing {key}: {e}")
            else:
                # try to serialize the value as a GeneralData
                try:
                    new_inputs[key] = GeneralData(value)
                except Exception as e:
                    raise ValueError(f"Error in serializing {key}: {e}")

    return new_inputs
