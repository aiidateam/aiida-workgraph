from typing import List, Union, Any
from aiida_workgraph.property import TaskProperty
from aiida import orm


class PropertyAny(TaskProperty):
    """A new class for Any type."""

    identifier: str = "workgraph.any"

    def validate(self, _: any) -> None:
        """No validation needed."""


class PropertyInt(TaskProperty):
    """A new class for integer type."""

    identifier: str = "workgraph.int"
    allowed_types = (int, str, type(None))


class PropertyFloat(TaskProperty):
    """A new class for float type."""

    identifier: str = "workgraph.float"
    allowed_types = (int, float, str, type(None))


class PropertyBool(TaskProperty):
    """A new class for bool type."""

    identifier: str = "workgraph.bool"
    allowed_types = (bool, int, str, type(None))


class PropertyString(TaskProperty):
    """A new class for string type."""

    identifier: str = "workgraph.string"
    allowed_types = (str, type(None))

    def validate(self, value: any) -> None:
        if not isinstance(value, self.allowed_types):
            raise TypeError(
                f"Expected value of type {self.allowed_types}, got {type(value).__name__} instead."
            )


class PropertyAiiDAInt(TaskProperty):
    """A new class for integer type."""

    identifier: str = "workgraph.aiida_int"
    allowed_types = (int, str, orm.Int, type(None))

    def set_value(self, value: Union[int, orm.Int, str] = None) -> None:
        if isinstance(value, int):
            value = orm.Int(value)
        super().set_value(value)


class PropertyAiiDAFloat(TaskProperty):
    """A new class for float type."""

    identifier: str = "workgraph.aiida_float"
    allowed_types = (int, float, orm.Int, orm.Float, str, type(None))

    def set_value(
        self, value: Union[int, float, orm.Int, orm.Float, str] = None
    ) -> None:
        if isinstance(value, (int, float)):
            value = orm.Float(value)
        super().set_value(value)


class PropertyAiiDABool(TaskProperty):
    """A new class for bool type."""

    identifier: str = "workgraph.aiida_bool"
    allowed_types = (int, bool, orm.Int, orm.Bool, str, type(None))

    def set_value(self, value: Union[int, bool, orm.Int, orm.Bool, str] = None) -> None:
        if isinstance(value, int):
            value = orm.Bool(value)
        super().set_value(value)


class PropertyAiiDAString(TaskProperty):
    """A new class for string type."""

    identifier: str = "workgraph.aiida_string"
    allowed_types = (int, float, str, orm.Str, type(None))

    def set_value(self, value: Union[int, float, str, orm.Str] = None) -> None:
        if isinstance(value, (int, float)):
            value = orm.Str(value)
        super().set_value(value)

    def validate(self, value: any) -> None:
        if not isinstance(value, self.allowed_types):
            raise TypeError(
                f"Expected value of type {self.allowed_types}, got {type(value).__name__} instead."
            )


class PropertyAiiDADict(TaskProperty):
    """A new class for Dict type."""

    identifier: str = "workgraph.aiida_dict"
    allowed_types = (dict, orm.Dict, str, type(None))

    def set_value(self, value: Union[dict, orm.Dict, str] = None) -> None:
        if isinstance(value, (dict)):
            value = orm.Dict(value)
        super().set_value(value)


# ====================================
class PropertyVector(TaskProperty):
    """Vector property"""

    identifier: str = "workgraph.vector"
    allowed_item_types = (object, type(None))

    def __init__(self, name, description="", size=3, default=None, update=None) -> None:
        self.size = size
        default = [] if default is None else default
        super().__init__(name, description, default, update)

    def validate(self, value: Any) -> None:
        """Validate the given value based on allowed types."""
        if value is not None:
            if len(value) != self.size:
                raise ValueError(
                    f"Invalid size: Expected {self.size}, got {len(value)} instead."
                )
            for item in value:
                if not isinstance(item, self.allowed_item_types):
                    raise ValueError(
                        f"Invalid item type: Expected {self.allowed_item_types}, got {type(item)} instead."
                    )

        super().validate(value)

    def set_value(self, value: List) -> None:
        self.validate(value)
        self._value = value
        if self.update:
            self.update()

    def copy(self):
        p = self.__class__(
            self.name, self.description, self.size, self.value, self.update
        )
        p.value = self.value
        return p

    def get_metadata(self):
        metadata = {"default": self.default, "size": self.size}
        return metadata


class PropertyAiiDAIntVector(PropertyVector):
    """A new class for integer vector type."""

    identifier: str = "workgraph.aiida_int_vector"
    allowed_types = (list, orm.List, type(None))
    allowed_item_types = (int, str, type(None))


class PropertyAiiDAFloatVector(PropertyVector):
    """A new class for float vector type."""

    identifier: str = "workgraph.aiida_float_vector"
    allowed_types = (list, orm.List, type(None))
    allowed_item_types = (int, float, str, type(None))


class PropertyStructureData(TaskProperty):
    """A new class for Any type."""

    identifier: str = "workgraph.aiida_structuredata"
    allowed_types = (orm.StructureData, str, type(None))


__all__ = [
    PropertyAny,
    PropertyAiiDAInt,
    PropertyAiiDAFloat,
    PropertyAiiDABool,
    PropertyAiiDAString,
    PropertyAiiDADict,
    PropertyAiiDAIntVector,
    PropertyAiiDAFloatVector,
]
