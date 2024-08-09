from typing import Dict, List, Union, Callable
from node_graph.property import NodeProperty
from node_graph.serializer import SerializeJson
from node_graph.properties.builtins import PropertyVector, PropertyAny
from aiida import orm


class PropertyAiiDAInt(NodeProperty, SerializeJson):
    """A new class for integer type."""

    identifier: str = "workgraph.aiida_int"
    data_type = "Int"

    def __init__(
        self,
        name: str,
        description: str = "",
        default: Union[int, None] = None,
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value: Union[int, orm.Int, str]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif isinstance(value, int):
            self._value = orm.Int(value)
            if self.update is not None:
                self.update()
        elif isinstance(value, orm.Int):
            self._value = value
            if self.update is not None:
                self.update()
        elif (
            isinstance(value, str)
            and value.rstrip().startswith("{{")
            and value.endswith("}}")
        ):
            self._value = value
        else:
            raise Exception("{} is not a integer.".format(value))


class PropertyAiiDAFloat(NodeProperty, SerializeJson):
    """A new class for float type."""

    identifier: str = "workgraph.aiida_float"
    data_type = "Float"

    def __init__(
        self,
        name: str,
        description: str = "",
        default: Union[int, None] = None,
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value: Union[float, orm.Float, int, orm.Int, str]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif isinstance(value, (int, float)):
            self._value = orm.Float(value)
            if self.update is not None:
                self.update()
        elif isinstance(value, (orm.Int, orm.Float)):
            self._value = value
            if self.update is not None:
                self.update()
        elif (
            isinstance(value, str)
            and value.rstrip().startswith("{{")
            and value.endswith("}}")
        ):
            self._value = value
        else:
            raise Exception("{} is not a float.".format(value))


class PropertyAiiDABool(NodeProperty, SerializeJson):
    """A new class for bool type."""

    identifier: str = "workgraph.aiida_bool"
    data_type = "Bool"

    def __init__(
        self,
        name: str,
        description: str = "",
        default: Union[bool, None] = None,
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value: Union[bool, orm.Bool, int, str]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif isinstance(value, (bool, int)):
            self._value = orm.Bool(value)
            if self.update is not None:
                self.update()
        elif isinstance(value, (orm.Bool)):
            self._value = value
            if self.update is not None:
                self.update()
        elif (
            isinstance(value, str)
            and value.rstrip().startswith("{{")
            and value.endswith("}}")
        ):
            self._value = value
        else:
            raise Exception("{} is not a bool.".format(value))


class PropertyAiiDAString(NodeProperty, SerializeJson):
    """A new class for string type."""

    identifier: str = "workgraph.aiida_string"
    data_type = "String"

    def __init__(
        self,
        name: str,
        description: str = "",
        default: Union[str, orm.Str, None] = None,
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value: Union[str, orm.Str]) -> None:
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif isinstance(value, str):
            self._value = orm.Str(value)
            if self.update is not None:
                self.update()
        elif isinstance(value, orm.Str):
            self._value = value
            if self.update is not None:
                self.update()
        elif (
            isinstance(value, str)
            and value.rstrip().startswith("{{")
            and value.endswith("}}")
        ):
            self._value = value
        else:
            raise Exception("{} is not a string.".format(value))


class PropertyAiiDADict(NodeProperty, SerializeJson):
    """A new class for Dict type."""

    identifier: str = "workgraph.aiida_dict"
    data_type = "Dict"

    def __init__(
        self,
        name: str,
        description: str = "",
        default: Union[dict, orm.Dict, None] = True,
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value: Union[Dict, orm.Dict, str]) -> None:
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif isinstance(value, dict):
            self._value = orm.Dict(value)
            if self.update is not None:
                self.update()
        elif isinstance(value, orm.Dict):
            self._value = value
            if self.update is not None:
                self.update()
        elif (
            isinstance(value, str)
            and value.rstrip().startswith("{{")
            and value.endswith("}}")
        ):
            self._value = value
        else:
            raise Exception("{} is not a dict.".format(value))


class PropertyAiiDAIntVector(PropertyVector):
    """A new class for integer vector type."""

    identifier: str = "workgraph.aiida_int_vector"
    data_type = "AiiDAIntVector"

    def __init__(
        self,
        name: str,
        description: str = "",
        size: int = 3,
        default: List[int] = [0, 0, 0],
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value: List[int]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif len(value) == self.size:
            for i in range(self.size):
                if isinstance(value[i], int):
                    self._value[i] = value[i]
                    if self.update is not None:
                        self.update()
                else:
                    raise Exception(
                        f"Set property {self.name} failed. {value[i]} is not a integer."
                    )
        else:
            raise Exception(
                "Length {} is not equal to the size {}.".format(len(value), self.size)
            )


class PropertyAiiDAFloatVector(PropertyVector):
    """A new class for float vector type."""

    identifier: str = "workgraph.aiida_float_vector"
    data_type = "AiiDAFloatVector"

    def __init__(
        self,
        name: str,
        description: str = "",
        size: int = 3,
        default: List[float] = [0.0, 0.0, 0.0],
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value: List[Union[int, float]]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif len(value) == self.size:
            for i in range(self.size):
                if isinstance(value[i], (int, float)):
                    self._value[i] = value[i]
                    if self.update is not None:
                        self.update()
                else:
                    raise Exception("{} is not a float.".format(value[i]))
        else:
            raise Exception(
                "Length {} is not equal to the size {}.".format(len(value), self.size)
            )

    def get_metadata(self):
        metadata = {"default": self.default, "size": self.size}
        return metadata


# ====================================
# Vector


class PropertyBoolVector(PropertyVector):
    """A new class for bool vector type."""

    identifier: str = "workgraph.bool_vector"
    data_type = "BoolVector"

    def __init__(
        self,
        name: str,
        description: str = "",
        size: int = 3,
        default: List[bool] = [False, False, False],
        update: Callable = None,
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value: List[Union[bool, int]]) -> None:
        # run the callback function
        if value is None:
            self._value = value
            if self.update is not None:
                self.update()
        elif len(value) == self.size:
            for i in range(self.size):
                if isinstance(value[i], (bool, int)):
                    self._value[i] = value[i]
                    if self.update is not None:
                        self.update()
                else:
                    raise Exception("{} is not a bool.".format(value[i]))
        else:
            raise Exception(
                "Length {} is not equal to the size {}.".format(len(value), self.size)
            )


__all__ = [
    PropertyAny,
    PropertyAiiDAInt,
    PropertyAiiDAFloat,
    PropertyAiiDABool,
    PropertyAiiDAString,
    PropertyAiiDADict,
    PropertyAiiDAIntVector,
    PropertyAiiDAFloatVector,
    PropertyBoolVector,
]
