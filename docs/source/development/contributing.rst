======================
Contributing
======================

Pre-commit
==========

To contribute to this repository, please enable pre-commit to ensure that the
code in commits conforms to the standards.

.. code-block:: console

    $ pip install -e .[tests, pre-commit]
    $ pre-commit install


Running tests
=============

To run all tests, use:

.. code-block:: console

    pytest


Setting path for Python executable for pythonjob tests
------------------------------------------------------

By default the pythonjob will use the executable ``python3`` to execute the
calcjobs in the tests. If you want to specify a different Python path
(e.g. from your environment manager), you can set the environment variable:

.. code-block:: console

    PYTEST_PYTHONJOB_PYTHON_EXEC_PATH=/home/user/pyvenv/workgraph-dev/bin/python pytest tests/test_python.py


Building the docs
=================

We use Sphinx to build the documentation. You need the requirements in the
extra ``.[docs]`` dependency and the ``docs/requirements.txt`` file.
We have a ``docs/Makefile`` that runs sphinx-build to build the docs.

.. code-block:: console

    pip install .[docs]
    pip install -r docs/requirements.txt
    make -C docs html
    <YOUR-BROWSER> docs/build/html/index.html


Creating a new Sphinx source file with executable code
------------------------------------------------------

We use **sphinx-gallery** to integrate executable code into the documentation.
For that we create a sphinx-gallery script (an extended Python file that can be
parsed by sphinx-gallery to generate an ``.rst`` file with more structure)
instead of a plain ``.rst`` file.

You can create a sphinx-gallery script from a Jupyter notebook using IDE like VSCode.

We put the converted sphinx-gallery script file ``<SCRIPT>.py`` into the gallery
source folder ``docs/gallery/<FOLDER>/autogen``, where ``<FOLDER>`` is the
folder you want to attach the generated file to in the sphinx source
(``docs/source``).

Then in the ``docs/source/<FOLDER>/index.rst`` add ``autogen/<SCRIPT>`` to the
toctree. The sphinx-gallery script will be converted to an ``.rst`` file during
build time.

To modify existing documentation, change the corresponding Python files under
``docs/gallery``.
