from typing import Any, Type
from node_graph.socket import NodeSocket
from aiida_workgraph.property import TaskProperty


class TaskSocket(NodeSocket):
    """Represent a socket of a Task in the AiiDA WorkGraph."""

    # use TaskProperty from aiida_workgraph.property
    # to override the default NodeProperty from node_graph
    node_property = TaskProperty

    @property
    def node_value(self):
        if not hasattr(self, '_node_value'):
            self._node_value = self.get_node_value()
        return self._node_value

    def get_node_value(self):
        "Directly return or set the _node_value attribute."

        # _node_value was set before already
        if hasattr(self, '_node_value'):
            pass

        # If data associated with Socket is AiiDA ORM, return again its value
        # Check for the nested case before, otherwise Socket value is matched
        # first and the ORM instance returned
        elif hasattr(self.value, "value"):
            # TODO: One could also check for isinstance of AiiDA ORM, however, that adds another import, and not sure if
            # TODO: it's really necessary here

            self._node_value = self.value.value
            # TODO: Possibly check here for AttributeDict, for which we return the `get_dict` result directly
            # TODO: However, with a specific `SocketAiiDADict` class, we could just overwrite the method there

        # If not, e.g. when Python base types are used, directly return the variable's value
        elif hasattr(self, "value"):
            self._node_value = self.value

        # TODO: This shouldn't even be a case? Shouldn't a `Socket` _always_ have an associated value, even if it's
        # TODO: `None` on instantiation
        else:
            msg = f"Socket <{self}> does not have an asociated `value` attribute."
            raise AttributeError(msg)

        return self._node_value

def build_socket_from_AiiDA(DataClass: Type[Any]) -> Type[TaskSocket]:
    """Create a socket class from AiiDA DataClass."""

    class AiiDATaskSocket(TaskSocket):
        """AiiDA Task Socket."""

        identifier: str = DataClass.__name__

        def __init__(
            self,
            name: str,
            parent: Any = None,
            type: str = "INPUT",
            index: int = 0,
            uuid: str = None,
            **kwargs: Any
        ) -> None:
            super().__init__(name, parent, type, index, uuid=uuid)
            self.add_property(DataClass, name, **kwargs)

        def get_serialize(self) -> dict:
            serialize = {"path": "aiida.orm.utils.serialize", "name": "serialize"}
            return serialize

        def get_deserialize(self) -> dict:
            deserialize = {
                "path": "aiida.orm.utils.serialize",
                "name": "deserialize_unsafe",
            }
            return deserialize

    return AiiDATaskSocket
