"""`Data` sub class to represent any data using pickle."""

from aiida import orm
import sys
import cloudpickle


class Dict(orm.Dict):
    @property
    def value(self):
        return self.get_dict()


class List(orm.List):
    @property
    def value(self):
        return self.get_list()


class GeneralData(orm.Data):
    """Data to represent a pickled value using cloudpickle."""

    FILENAME = "value.pkl"  # Class attribute to store the filename

    def __init__(self, value=None, **kwargs):
        """Initialize a `GeneralData` node instance.

        :param value: raw Python value to initialize the `GeneralData` node from.
        """
        super().__init__(**kwargs)
        self.set_value(value)

    def __str__(self):
        return f"{super().__str__()} : {self.get_value()}"

    @property
    def value(self):
        """Return the contents of this node.

        :return: The unpickled value.
        """
        return self.get_value()

    @value.setter
    def value(self, value):
        self.set_value(value)

    def get_value(self):
        """Return the contents of this node, unpickling the stored value.

        :return: The unpickled value.
        """
        self._check_pickle_protocol()
        return self._get_value_from_file()

    def _get_value_from_file(self):
        """Read the pickled value from file and return it."""
        with self.base.repository.open(self.FILENAME, mode="rb") as f:
            return cloudpickle.loads(f.read())  # Deserialize the value

    def set_value(self, value):
        """Set the contents of this node by pickling the provided value.

        :param value: The Python value to pickle and store.
        """
        # Serialize the value and store it
        serialized_value = cloudpickle.dumps(value)
        self.base.repository.put_object_from_bytes(serialized_value, self.FILENAME)

        # Store relevant metadata
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        self.base.attributes.set("python_version", python_version)
        self.base.attributes.set("serializer_module", cloudpickle.__name__)
        self.base.attributes.set("serializer_version", cloudpickle.__version__)
        self.base.attributes.set("pickle_protocol", cloudpickle.DEFAULT_PROTOCOL)

    def _check_pickle_protocol(self):
        """Check if the stored pickle protocol matches the current environment.

        Raises:
            ImportError: If the pickle protocol does not match.
        """
        stored_pickle_protocol = self.base.attributes.get("pickle_protocol")

        if stored_pickle_protocol != cloudpickle.DEFAULT_PROTOCOL:
            raise ImportError(
                f"Incompatible pickle protocol: expected {stored_pickle_protocol}, "
                f"but found {cloudpickle.DEFAULT_PROTOCOL}. Please ensure the correct protocol is used."
            )
