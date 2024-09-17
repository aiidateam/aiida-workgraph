from collections.abc import Iterator 

from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.nodes.data.data import Data
import copy
from aiida.orm.nodes.data import to_aiida_type

__all__ = ('Iterator',)


class AiidaIterator(Data, Iterator):
    """`Data` sub class to represent a iterator."""

    _ITERATOR_KEY = 'iter'

    def __init__(self, value=None, **kwargs):
        """Initialise a ``Iterator`` node instance.

        :param value: iterator to initialise the ``Iterator`` node from
        """
        data = value or kwargs.pop('iter', [])
        super().__init__(**kwargs)
        self.set_iterator(data)

    def __next__(self):
        iterator = self.get_iterator()
        iterator.__next__()
        if not self._using_reference():
            self.set_iterator(iterator)

    def get_iterator(self):
        """Return the contents of this node.

        :return: a iterator
        """
        return self.base.attributes.get(self._ITERATOR_KEY)

    def set_iterator(self, data):
        """Set the contents of this node.

        :param data: the iterator to set
        """
        if not hasattr(data, "__iter__"):
            raise TypeError('Must supply type that implements __iter__')
        self.base.attributes.set(self._ITERATOR_KEY, copy.deepcopy(data))

    def _using_reference(self):
        """This function tells the class if we are using a iterator reference.  This
        means that calls to self.get_iterator return a reference rather than a copy
        of the underlying iterator and therefore self.set_iterator need not be called.
        This knwoledge is essential to make sure this class is performant.

        Currently the implementation assumes that if the node needs to be
        stored then it is using the attributes cache which is a reference.

        :return: True if using self.get_iterator returns a reference to the
            underlying sequence.  False otherwise.
        :rtype: bool
        """
        return not self.is_stored
        
from _collections_abc import list_iterator
@to_aiida_type.register(list_iterator)
def _(value):
    return AiidaIterator(value)
