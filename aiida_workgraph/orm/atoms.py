# -*- coding: utf-8 -*-
###########################################################################
# Copyright (c), The AiiDA team. All rights reserved.                     #
# This file is part of the AiiDA code.                                    #
#                                                                         #
# The code is hosted on GitHub at https://github.com/aiidateam/aiida-core #
# For further information on the license, see the LICENSE.txt file        #
# For further information please visit http://www.aiida.net               #
###########################################################################
"""`Data` sub class to represent a list."""

from aiida.orm import Data
from ase import Atoms

__all__ = ("AtomsData",)


class AtomsData(Data):
    """`Data to represent a ASE Atoms."""

    _cached_atoms = None

    def __init__(self, value=None, **kwargs):
        """Initialise a ``List`` node instance.

        :param value: list to initialise the ``List`` node from
        """
        data = value or kwargs.pop("atoms", Atoms())
        super().__init__(**kwargs)
        self.set_atoms(data)

    @property
    def value(self):
        return self.get_atoms()

    def initialize(self):
        super().initialize()
        self._cached_atoms = None

    def __getitem__(self, item):
        return self.get_atoms()[item]

    def __setitem__(self, key, value):
        data = self.get_atoms()
        data[key] = value
        if not self._using_atoms_reference():
            self.set_atoms(data)

    def __delitem__(self, key):
        data = self.get_atoms()
        del data[key]
        if not self._using_atoms_reference():
            self.set_atoms(data)

    def __len__(self):
        return len(self.get_atoms())

    def __str__(self):
        return f"{super().__str__()} : {self.get_atoms()}"

    def __eq__(self, other):
        if isinstance(other, Atoms):
            return self.get_atoms() == other.get_atoms()
        return self.get_atoms() == other

    def append(self, value):
        data = self.get_atoms()
        data.append(value)
        if not self._using_atoms_reference():
            self.set_atoms(data)

    def extend(self, value):  # pylint: disable=arguments-renamed
        data = self.get_atoms()
        data.extend(value)
        if not self._using_atoms_reference():
            self.set_atoms(data)

    def get_atoms(self):
        """Return the contents of this node.

        :return: a Atoms
        """
        import pickle

        def get_atoms_from_file(self):
            filename = "atoms.pkl"
            # Open a handle in binary read mode as the arrays are written as binary files as well
            with self.base.repository.open(filename, mode="rb") as f:
                return pickle.loads(f.read())  # pylint: disable=unexpected-keyword-arg

        # Return with proper caching if the node is stored, otherwise always re-read from disk
        if not self.is_stored:
            return get_atoms_from_file(self)

        if self._cached_atoms is None:
            self._cached_atoms = get_atoms_from_file(self)

        return self._cached_atoms

    def set_atoms(self, atoms):
        """Set the contents of this node.

        :param atoms: the atoms to set
        """
        import pickle

        if not isinstance(atoms, Atoms):
            raise TypeError("Must supply Atoms type")
        self.base.repository.put_object_from_bytes(pickle.dumps(atoms), "atoms.pkl")
        formula = atoms.get_chemical_formula()
        # Store the array name and shape for querying purposes
        self.base.attributes.set("formula", formula)

    def _using_atoms_reference(self):
        """
        This function tells the class if we are using a list reference.  This
        means that calls to self.get_atoms return a reference rather than a copy
        of the underlying list and therefore self.set_atoms need not be called.
        This knwoledge is essential to make sure this class is performant.

        Currently the implementation assumes that if the node needs to be
        stored then it is using the attributes cache which is a reference.

        :return: True if using self.get_atoms returns a reference to the
            underlying sequence.  False otherwise.
        :rtype: bool
        """
        return self.is_stored
