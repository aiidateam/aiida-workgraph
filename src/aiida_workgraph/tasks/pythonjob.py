from typing import Any, Dict
from aiida import orm
from aiida_pythonjob.data.serializer import general_serializer
from aiida_workgraph.task import Task


class PythonJob(Task):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def update_from_dict(self, data: Dict[str, Any], **kwargs) -> "PythonJob":
        """Overwrite the update_from_dict method to handle the PythonJob data."""
        self.deserialize_pythonjob_data(data["inputs"])
        super().update_from_dict(data)

    @classmethod
    def serialize_pythonjob_data(
        cls, input_data: Dict[str, Any], is_function_input: bool = False
    ) -> None:
        """Serialize the properties for PythonJob."""

        for input in input_data.values():
            if is_function_input or input["metadata"].get("is_function_input", False):
                if input["identifier"] == "workgraph.namespace":
                    cls.serialize_pythonjob_data(
                        input["sockets"], is_function_input=True
                    )
                elif input.get("property", {}).get("value") is not None:
                    cls.serialize_socket_data(input)

    @classmethod
    def deserialize_pythonjob_data(
        cls, input_data: Dict[str, Any], is_function_input: bool = False
    ) -> None:
        """
        Process the task data dictionary for a PythonJob.
        It load the orignal Python data from the AiiDA Data node for the
        args and kwargs of the function.

        Args:
            tdata (Dict[str, Any]): The input data dictionary.

        Returns:
            Dict[str, Any]: The processed data dictionary.
        """

        for input in input_data.values():
            if is_function_input or input["metadata"].get("is_function_input", False):
                if input["identifier"] == "workgraph.namespace":
                    print("deserialize namespace: ", input["name"])
                    cls.deserialize_pythonjob_data(
                        input["sockets"], is_function_input=True
                    )
                else:
                    print("deserialize socket: ", input["name"])
                    cls.deserialize_socket_data(input)

    @classmethod
    def serialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        value = data.get("property", {}).get("value")
        if value is None or isinstance(value, orm.Data):
            return
        data["property"]["value"] = general_serializer(value)

    @classmethod
    def deserialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        value = data.get("property", {}).get("value")
        if isinstance(value, orm.Data):
            data["property"]["value"] = value.value
