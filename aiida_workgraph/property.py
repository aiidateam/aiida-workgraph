from node_graph.property import NodeProperty as GraphNodeProperty


class NodeProperty(GraphNodeProperty):
    """Represent a property of a Node in the AiiDA WorkGraph."""

    @classmethod
    def new(cls, identifier, name=None, data={}):
        """Create a property from a identifier."""
        # use property_pool from aiida_workgraph.properties
        # to override the default property_pool from node_graph
        from aiida_workgraph.properties import property_pool

        # build the node on the fly if the identifier is a callable
        if callable(identifier):
            identifier = build_property_from_AiiDA(identifier)
        return super().new(
            identifier, name=name, data=data, property_pool=property_pool
        )


def build_property_from_AiiDA(DataClass):
    """Create a property class from AiiDA DataClass."""

    class AiiDANodeProperty(NodeProperty):
        identifier: str = DataClass.__name__

        def set_value(self, value):
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

        def get_serialize(self):
            serialize = {"path": "aiida.orm.utils.serialize", "name": "serialize"}
            return serialize

        def get_deserialize(self):
            deserialize = {
                "path": "aiida.orm.utils.serialize",
                "name": "deserialize_unsafe",
            }
            return deserialize

    return AiiDANodeProperty
