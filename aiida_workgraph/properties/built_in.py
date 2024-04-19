from node_graph.property import NodeProperty
from node_graph.serializer import SerializeJson, SerializePickle
from node_graph.properties.builtin import (
    VectorProperty,
    BaseDictProperty,
    BaseListProperty,
    IntProperty,
    BoolProperty,
    FloatProperty,
    StringProperty,
)
from aiida import orm


class GeneralProperty(NodeProperty, SerializePickle):
    """A new class for General type."""

    identifier: str = "General"
    data_type = "General"

    def __init__(self, name, description="", default=None, update=None) -> None:
        super().__init__(name, description, default, update)


class AiiDAIntProperty(NodeProperty, SerializeJson):
    """A new class for integer type."""

    identifier: str = "AiiDAInt"
    data_type = "Int"

    def __init__(self, name, description="", default=0, update=None) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value):
        # run the callback function
        if isinstance(value, int):
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


class AiiDAFloatProperty(NodeProperty, SerializeJson):
    """A new class for float type."""

    identifier: str = "AiiDAFloat"
    data_type = "Float"

    def __init__(self, name, description="", default=0.0, update=None) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value):
        # run the callback function
        if isinstance(value, (int, float)):
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


class AiiDABoolProperty(NodeProperty, SerializeJson):
    """A new class for bool type."""

    identifier: str = "AiiDABool"
    data_type = "Bool"

    def __init__(self, name, description="", default=True, update=None) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value):
        # run the callback function
        if isinstance(value, (bool, int)):
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


class AiiDAStringProperty(NodeProperty, SerializeJson):
    """A new class for string type."""

    identifier: str = "AiiDAString"
    data_type = "String"

    def __init__(self, name, description="", default=None, update=None) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value):
        if isinstance(value, str):
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


class AiiDADictProperty(NodeProperty, SerializeJson):
    """A new class for Dict type."""

    identifier: str = "AiiDADict"
    data_type = "Dict"

    def __init__(self, name, description="", default=None, update=None) -> None:
        super().__init__(name, description, default, update)

    def set_value(self, value):
        if isinstance(value, dict):
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


class AiiDAIntVectorProperty(VectorProperty):
    """A new class for integer vector type."""

    identifier: str = "AiiDAIntVector"
    data_type = "AiiDAIntVector"

    def __init__(
        self, name, description="", size=3, default=[0, 0, 0], update=None
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value):
        # run the callback function
        if len(value) == self.size:
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


class AiiDAFloatVectorProperty(VectorProperty):
    """A new class for float vector type."""

    identifier: str = "AiiDAFloatVector"
    data_type = "AiiDAFloatVector"

    def __init__(
        self, name, description="", size=3, default=[0, 0, 0], update=None
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value):
        # run the callback function
        if len(value) == self.size:
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


class BoolVectorProperty(VectorProperty):
    """A new class for bool vector type."""

    identifier: str = "BoolVector"
    data_type = "BoolVector"

    def __init__(
        self, name, description="", size=3, default=[0, 0, 0], update=None
    ) -> None:
        super().__init__(name, description, size, default, update)

    def set_value(self, value):
        # run the callback function
        if len(value) == self.size:
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


property_list = [
    IntProperty,
    FloatProperty,
    BoolProperty,
    StringProperty,
    GeneralProperty,
    BaseDictProperty,
    BaseListProperty,
    AiiDAIntProperty,
    AiiDAFloatProperty,
    AiiDAStringProperty,
    AiiDABoolProperty,
    AiiDADictProperty,
    AiiDAIntVectorProperty,
    AiiDAFloatVectorProperty,
]
