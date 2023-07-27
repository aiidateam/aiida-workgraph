============
Installation
============

.. _installation:requirements:

Requirements
============

To work with ``aiida-worktree``, you should have:

* installed ``aiida-core``
* configured an AiiDA profile.

Please refer to the `documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html>`_ of ``aiida-core`` for detailed instructions.


.. _installation:installation:

Installation
============


The recommended method of installation is to use the Python package manager |pip|_:

.. code-block:: console

    $ pip install aiida-worktree

This will install the latest stable version that was released to PyPI.

To install the package from source, first clone the repository and then install using |pip|_:

.. code-block:: console

    $ git clone https://github.com/superstar54/aiida-worktree
    $ pip install -e aiida-worktree

The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.



.. |pip| replace:: ``pip``
.. _pip: https://pip.pypa.io/en/stable/
