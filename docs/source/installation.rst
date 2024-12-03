============
Installation
============

.. _installation:requirements:

Requirements
============

To work with ``aiida-workgraph``, you should have:

* installed ``aiida-core``
* configured an AiiDA profile.

Please refer to the `documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html>`_ of ``aiida-core`` for detailed instructions.


.. _installation:installation:

Installation
============


The recommended method of installation is to use the Python package manager |pip|_:

.. code-block:: console

    $ pip install aiida-workgraph
    $ # install web ui package if you want to use the web ui
    $ pip install aiida-workgraph-web-ui

This will install the latest stable version that was released to PyPI.

To install the package from source, first clone the repository and then install using |pip|_:

.. code-block:: console

    $ git clone https://github.com/aiidateam/aiida-workgraph
    $ cd aiida-workgraph
    $ pip install -e .

The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.
To install the web app you need to in addition build the JavaScript packages:

.. code-block:: console

    $ cd aiida-workgraph
    $ pip install -e .


.. |pip| replace:: ``pip``
.. _pip: https://pip.pypa.io/en/stable/
