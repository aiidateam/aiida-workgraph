from typing import Any, Type, Union, Optional, Callable
from node_graph.property import NodeProperty


class TaskProperty(NodeProperty):
    """Represent a property of a Task in the AiiDA WorkGraph."""

    def validate(self, value: any) -> None:
        super().validate(value)

    @classmethod
    def new(
        cls, identifier: Union[Callable, str], name: Optional[str] = None, **kwargs
    ) -> "TaskProperty":
        """Create a property from a identifier."""
        # use PropertyPool from aiida_workgraph.properties
        # to override the default PropertyPool from node_graph
        from aiida_workgraph.properties import PropertyPool

        # build the task on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_property_from_AiiDA(identifier)
        return super().new(identifier, name=name, PropertyPool=PropertyPool, **kwargs)


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
