"""`Data` sub class to represent any data using pickle."""

from aiida import orm


class Dict(orm.Dict):
    @property
    def value(self):
        return self.get_dict()


class List(orm.List):
    @property
    def value(self):
        return self.get_list()


class GeneralData(orm.Data):
    """`Data to represent a pickled value."""

    def __init__(self, value=None, **kwargs):
        """Initialise a ``General`` node instance.

        :param value: list to initialise the ``List`` node from
        """
        super().__init__(**kwargs)
        self.set_value(value)

    def __str__(self):
        return f"{super().__str__()} : {self.get_value()}"

    @property
    def value(self):
        """Return the contents of this node.

        :return: a value
        """
        return self.get_value()

    @value.setter
    def value(self, value):
        return self.set_value(value)

    def get_value(self):
        """Return the contents of this node.

        :return: a value
        """
        import cloudpickle as pickle

        def get_value_from_file(self):
            filename = "value.pkl"
            # Open a handle in binary read mode as the arrays are written as binary files as well
            with self.base.repository.open(filename, mode="rb") as f:
                return pickle.loads(f.read())  # pylint: disable=unexpected-keyword-arg

        # Return with proper caching if the node is stored, otherwise always re-read from disk
        return get_value_from_file(self)

    def set_value(self, value):
        """Set the contents of this node.

        :param value: the value to set
        """
        import cloudpickle as pickle
        import sys

        self.base.repository.put_object_from_bytes(pickle.dumps(value), "value.pkl")
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        self.base.attributes.set("python_version", python_version)

    def _using_value_reference(self):
        """This function tells the class if we are using a list reference.  This
        means that calls to self.get_value return a reference rather than a copy
        of the underlying list and therefore self.set_value need not be called.
        This knwoledge is essential to make sure this class is performant.

        Currently the implementation assumes that if the node needs to be
        stored then it is using the attributes cache which is a reference.

        :return: True if using self.get_value returns a reference to the
            underlying sequence.  False otherwise.
        :rtype: bool
        """
        return self.is_stored
