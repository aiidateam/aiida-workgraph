from typing import Any, Dict, List
from aiida import orm
from aiida_workgraph.orm.serializer import general_serializer
from aiida_workgraph.task import Task


class PythonJob(Task):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    function_kwargs: List = None

    def update_from_dict(self, data: Dict[str, Any], **kwargs) -> "PythonJob":
        """Overwrite the update_from_dict method to handle the PythonJob data."""
        self.deserialize_pythonjob_data(data)
        self.function_kwargs = data.get("function_kwargs", [])
        super().update_from_dict(data)

    def to_dict(self) -> Dict[str, Any]:
        data = super().to_dict()
        data["function_kwargs"] = self.function_kwargs
        self.serialize_pythonjob_data(data)
        return data

    def serialize_pythonjob_data(self, tdata: Dict[str, Any]):
        """Serialize the properties for PythonJob."""

        input_kwargs = tdata.get("function_kwargs", [])
        for name in input_kwargs:
            tdata["inputs"][name]["property"]["value"] = self.serialize_socket_data(
                tdata["inputs"][name]
            )

    def deserialize_pythonjob_data(self, tdata: Dict[str, Any]) -> None:
        """
        Process the task data dictionary for a PythonJob.
        It load the orignal Python data from the AiiDA Data node for the
        args and kwargs of the function.

        Args:
            tdata (Dict[str, Any]): The input data dictionary.

        Returns:
            Dict[str, Any]: The processed data dictionary.
        """
        input_kwargs = tdata.get("function_kwargs", [])

        for name in input_kwargs:
            if name in tdata["inputs"]:
                tdata["inputs"][name]["property"][
                    "value"
                ] = self.deserialize_socket_data(tdata["inputs"][name])

    def find_input_socket(self, name):
        """Find the output with the given name."""
        if name in self.inputs:
            return self.inputs[name]
        return None

    def serialize_socket_data(self, data: Dict[str, Any]) -> Any:
        name = data["name"]
        if data.get("identifier", "Any").upper() == "WORKGRAPH.NAMESPACE":
            if isinstance(data["property"]["value"], dict):
                serialized_result = {}
                for key, value in data["property"]["value"].items():
                    full_name = f"{name}.{key}"
                    full_name_output = self.find_input_socket(full_name)
                    if (
                        full_name_output
                        and full_name_output.get("identifier", "Any").upper()
                        == "WORKGRAPH.NAMESPACE"
                    ):
                        serialized_result[key] = self.serialize_socket_data(
                            full_name_output
                        )
                    else:
                        serialized_result[key] = general_serializer(value)
                return serialized_result
            else:
                raise ValueError("Namespace socket should be a dictionary.")
        else:
            if isinstance(data["property"]["value"], orm.Data):
                return data["property"]["value"]
            return general_serializer(data["property"]["value"])

    def deserialize_socket_data(self, data: Dict[str, Any]) -> Any:
        name = data["name"]
        if data.get("identifier", "Any").upper() == "WORKGRAPH.NAMESPACE":
            if isinstance(data["property"]["value"], dict):
                deserialized_result = {}
                for key, value in data["property"]["value"].items():
                    full_name = f"{name}.{key}"
                    full_name_output = self.find_input_socket(full_name)
                    if (
                        full_name_output
                        and full_name_output.get("identifier", "Any").upper()
                        == "WORKGRAPH.NAMESPACE"
                    ):
                        deserialized_result[key] = self.deserialize_socket_data(
                            full_name_output
                        )
                    else:
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
