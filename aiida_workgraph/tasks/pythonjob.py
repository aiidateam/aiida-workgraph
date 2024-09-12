from typing import Any, Dict
from aiida import orm
from aiida_workgraph.orm.serializer import general_serializer
from aiida_workgraph.task import Task


class PythonJob(Task):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    @classmethod
    def get_function_kwargs(cls, data) -> Dict[str, Any]:
        input_kwargs = set()
        for name in data["kwargs"]:
            # all the kwargs are after computer is the input for the PythonJob, should be AiiDA Data node
            if name == "computer":
                break
            input_kwargs.add(name)
        return input_kwargs

    def update_from_dict(cls, data: Dict[str, Any], **kwargs) -> "PythonJob":
        """Overwrite the update_from_dict method to handle the PythonJob data."""
        cls.deserialize_pythonjob_data(data)
        return super().update_from_dict(data)

    def to_dict(self) -> Dict[str, Any]:
        data = super().to_dict()
        self.serialize_pythonjob_data(data)
        return data

    @classmethod
    def serialize_pythonjob_data(cls, tdata: Dict[str, Any]):
        """Serialize the properties for PythonJob."""

        input_kwargs = cls.get_function_kwargs(tdata)
        for name in input_kwargs:
            prop = tdata["properties"][name]
            # if value is not None, not {}
            if not (
                prop["value"] is None
                or isinstance(prop["value"], dict)
                and prop["value"] == {}
            ):
                prop["value"] = general_serializer(prop["value"])

    @classmethod
    def deserialize_pythonjob_data(cls, tdata: Dict[str, Any]) -> None:
        """
        Process the task data dictionary for a PythonJob.
        It load the orignal Python data from the AiiDA Data node for the
        args and kwargs of the function.

        Args:
            tdata (Dict[str, Any]): The input data dictionary.

        Returns:
            Dict[str, Any]: The processed data dictionary.
        """
        input_kwargs = cls.get_function_kwargs(tdata)

        for name in input_kwargs:
            if name in tdata["properties"]:
                value = tdata["properties"][name]["value"]
                if isinstance(value, orm.Data):
                    value = value.value
                elif value is not None and value != {}:
                    raise ValueError(f"There something wrong with the input {name}")
                tdata["properties"][name]["value"] = value
