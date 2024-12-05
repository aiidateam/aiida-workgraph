from typing import Any, Dict
from aiida import orm
from aiida_pythonjob.data.serializer import general_serializer
from aiida_workgraph.task import Task


class PythonJob(Task):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def update_from_dict(self, data: Dict[str, Any], **kwargs) -> "PythonJob":
        """Overwrite the update_from_dict method to handle the PythonJob data."""
        self.deserialize_pythonjob_data(data)
        super().update_from_dict(data)

    @classmethod
    def serialize_pythonjob_data(cls, tdata: Dict[str, Any]):
        """Serialize the properties for PythonJob."""

        for input in tdata["inputs"].values():
            if input["metadata"].get("is_function_input", False):
                input["property"]["value"] = cls.serialize_socket_data(input)

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

        for input in tdata["inputs"].values():
            if input["metadata"].get("is_function_input", False):
                input["property"]["value"] = cls.deserialize_socket_data(input)

    @classmethod
    def serialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        if data.get("identifier", "Any").upper() == "WORKGRAPH.NAMESPACE":
            if data["property"]["value"] is None:
                return None
            if isinstance(data["property"]["value"], dict):
                serialized_result = {}
                for key, value in data["property"]["value"].items():
                    serialized_result[key] = general_serializer(value)
                return serialized_result
            else:
                raise ValueError("Namespace socket should be a dictionary.")
        else:
            if isinstance(data["property"]["value"], orm.Data):
                return data["property"]["value"]
            return general_serializer(data["property"]["value"])

    @classmethod
    def deserialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        if data.get("identifier", "Any").upper() == "WORKGRAPH.NAMESPACE":
            if isinstance(data["property"]["value"], dict):
                deserialized_result = {}
                for key, value in data["property"]["value"].items():
                    if isinstance(value, orm.Data):
                        deserialized_result[key] = value.value
                    else:
                        deserialized_result[key] = value
                return deserialized_result
            else:
                raise ValueError("Namespace socket should be a dictionary.")
        else:
            if isinstance(data["property"]["value"], orm.Data):
                return data["property"]["value"].value
            return data["property"]["value"]
