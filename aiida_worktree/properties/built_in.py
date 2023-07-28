from scinode.core.property import NodeProperty
from scinode.serialization.built_in import SerializeJson, SerializePickle
from aiida import orm


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


property_list = [
    AiiDAIntProperty,
    AiiDAFloatProperty,
    AiiDAStringProperty,
    AiiDABoolProperty,
    AiiDADictProperty,
]
