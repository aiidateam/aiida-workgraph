from typing import Any, Type, Union, Dict, Optional, Callable
from node_graph.property import NodeProperty
import re


class TaskProperty(NodeProperty):
    """Represent a property of a Task in the AiiDA WorkGraph."""

    def validate(self, value: any) -> None:
        if isinstance(value, str) and not re.search(r"\{\{.*?\}\}", value):
            raise TypeError(
                f"""Expected value of type {self.allowed_types}, got {type(value).__name__} instead.
If you want to use variable from context, use double curly braces like this: {{{{variable_name}}}}"""
            )
        super().validate(value)

    @classmethod
    def new(
        cls,
        identifier: Union[Callable, str],
        name: Optional[str] = None,
        data: Dict[str, Any] = {},
    ) -> "TaskProperty":
        """Create a property from a identifier."""
        # use property_pool from aiida_workgraph.properties
        # to override the default property_pool from node_graph
        from aiida_workgraph.properties import property_pool

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_property_from_AiiDA(identifier)
        return super().new(
            identifier, name=name, data=data, property_pool=property_pool
        )


def build_property_from_AiiDA(DataClass: Type[Any]) -> Type[TaskProperty]:
    """Create a property class from AiiDA DataClass."""

    class AiiDATaskProperty(TaskProperty):
        identifier: str = DataClass.__name__

        def set_value(self, value: Any) -> None:
            # run the callback function
            if isinstance(value, DataClass):
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
                raise Exception("{} is not an {}.".format(value, DataClass.__name__))

    return AiiDATaskProperty
