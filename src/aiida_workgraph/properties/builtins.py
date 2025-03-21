from typing import List, Any
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
    allowed_types = (int, orm.Int, type(None))


class PropertyFloat(TaskProperty):
    """A new class for float type."""

    identifier: str = "workgraph.float"
    allowed_types = (int, float, orm.Int, orm.Float, type(None))


class PropertyBool(TaskProperty):
    """A new class for bool type."""

    identifier: str = "workgraph.bool"
    allowed_types = (bool, int, orm.Bool, orm.Int, type(None))


class PropertyString(TaskProperty):
    """A new class for string type."""

    identifier: str = "workgraph.string"
    allowed_types = (str, orm.Str, type(None))


class PropertyList(TaskProperty):
    """A new class for List type."""

    identifier: str = "workgraph.list"
    allowed_types = (list, orm.List, type(None))


class PropertyDict(TaskProperty):
    """A new class for Dict type."""

    identifier: str = "workgraph.dict"
    allowed_types = (dict, orm.Dict, type(None))


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
    allowed_item_types = (int, type(None))


class PropertyAiiDAFloatVector(PropertyVector):
    """A new class for float vector type."""

    identifier: str = "workgraph.aiida_float_vector"
    allowed_types = (list, orm.List, type(None))
    allowed_item_types = (int, float, type(None))


class PropertyStructureData(TaskProperty):
    """A new class for Any type."""

    identifier: str = "workgraph.aiida_structuredata"
    allowed_types = (orm.StructureData, type(None))
