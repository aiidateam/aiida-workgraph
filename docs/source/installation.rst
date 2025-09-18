============
Installation
============

.. _installation:requirements:

Requirements
____________

To work with ``aiida-workgraph``, you'll need:

* ``aiida-core`` (can also be installed directly with ``aiida-workgraph`` as it is a dependency)
* a working AiiDA profile (e.g., using AiiDA's ``verdi presto``)

Please refer to the `documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html>`_ of ``aiida-core`` for detailed instructions on AiiDA profile creation.


.. _installation:installation:

Using pip
_________

The recommended method of installation is to use the Python package manager |pip|_:

.. code-block:: console

    $ pip install aiida-workgraph

This will install the latest stable version that was released to `PyPI <https://pypi.org/project/aiida-workgraph>`_ (including all required dependencies).

From source
___________

To install the package from source, first clone the `repository <https://github.com/aiidateam/aiida-workgraph>`_ from GitHub and then install using |pip|_:

.. code-block:: console

    $ git clone https://github.com/aiidateam/aiida-workgraph
    $ cd aiida-workgraph
    $ pip install -e .

.. note::
   The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.
   Without this flag, the package is installed as a regular copy, and you would need to reinstall it each time you make changes to the source code.

.. _installation:gui:

GUI
___

In addition to the main ``aiida-workgraph`` Python package, we also provide a graphical user interface (GUI), which can be installed as follows:

.. code-block:: console

    $ pip install aiida-gui-workgraph

.. warning::
   **Experimental Feature**

   The GUI is still an experimental feature and under active development.
   Changes may be applied in future versions.
   Use with caution in production environments.

.. |pip| replace:: ``pip``
.. _pip: https://pip.pypa.io/en/stable/
