"""
Quick Start
===========

"""

# %%
# Introduction
# ------------
#
# This quick-start tutorial is intended to demonstrate the very basics of ``WorkGraph``.
# We do so through a simple arithmetic workflow.
# At the end of the tutorial, we will suggest some next steps to help you apply ``WorkGraph`` principles to your real-world workflows.

# %%
# Setup
# -----
#
# Let's first install ``aiida-workgraph`` (this will also install ``aiida-core`` automatically):
#
# .. code:: console
#
#    $ pip install aiida-workgraph
#
# To interact with the AiiDA database, we need to load an AiiDA profile.
# If you haven't configured one yet, you can do so by running the following command:
#
# .. code:: console
#
#    $ verdi presto
#
# .. note::
#
#    The rest of the quick start tutorial is best followed using a Jupyter notebook.
#
# To load your profile, run the following code:

from aiida import load_profile

load_profile()

# %%
# .. tip::
#
#    Some features of ``WorkGraph`` are best demonstrated interactively.
#    We recommend using the AiiDA GUI for this.
#    You can learn how to run and use the GUI in the :doc:`../gui/web` section.

# %%
# Simple workflow
# ---------------
#
# Suppose you want to compute ``(x + y) * z`` as a workflow of two distinct tasks: (i) addition, and (ii) multiplication.
# Let's see how we can do this using ``WorkGraph``:
#
# First, let's define our tasks using the ``@task`` decorator:

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# .. note::
#
#    To learn more about the ``@task`` decorator, visit the :doc:`../concept/autogen/task_concept` concept section.
#
# Next, we use our tasks to define a workflow using the ``@task.graph`` decorator:


@task.graph
def AddMultiply(x, y, z):
    the_sum = add(x=x, y=y).result
    return multiply(x=the_sum, y=z).result


# %%
#  .. note::
#
#     To learn more about the ``@task.graph`` decorator, visit the :doc:`../concept/autogen/graph_builder_concept` concept section.
#
# The workflow accepts three input parameters, ``x``, ``y``, and ``z``.
# It passes ``x`` and ``y`` to the ``add`` task, which computes their sum.
# The sum is then passed to the ``multiply`` task along with ``z``.
# The result of the multiplication task is returned as the output of the workflow.
#
# .. note::
#
#    Calling a ``@task``-decorated Python function returns a socket namespace - a collection of output sockets.
#    We do this to allow tasks to have multiple outputs.
#    We access a given socket with ``.<socket_name>``.
#    In the absence of explicit definition of output sockets, the default socket is named ``result``.
#
#    You can visit the following sections for more information about these features:
#
#    - :doc:`../concept/autogen/socket_concept` concept section
#    - :doc:`../howto/autogen/graph_level` how-to section
#
# Let's build the workflow with some input and visualize it:
#
# .. important::
#
#    If you are following along in a Jupyter notebook, omit the ``to_html()`` method call.

wg = AddMultiply.build_graph(x=2, y=3, z=4)
wg.to_html()

# %%
# We can see our two tasks, the assignment of the sum to the multiplication task, and the subsequent assignment of the product to the workflow (graph) result.
#
# Let's run the workflow with some inputs and inspect the result:

wg.run()

print("Result:", wg.outputs.result.value)

# %%
# We can also access the results of the individual tasks:

print("Result of addition:", wg.tasks.add.outputs.result.value)
print("Result of multiplication:", wg.tasks.multiply.outputs.result.value)

# %%
# .. note::
#
#    Functional tasks (i.e., Python functions decorated with ``@task``) are automatically named by their function name.
#    If you call the same task multiple times, subsequent calls will be suffixed with a number, e.g., ``add1``, ``add2``, etc.

# %%
# .. tip::
#
#    If you're running the AiiDA GUI, you can visualize the executed workflow interactively.
#    Click on the ``PK`` field of the submitted workflow (look for *WorkGraph<AddMultiply>*) to view its details.
#
# Let's have a look at the full provenance of our executed workflow:

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Summary
# -------
#
# This wraps up our quick start tutorial.
# You have learned how to:
#
# - Define tasks using the ``@task`` decorator
# - Define a workflow using the ``@task.graph`` decorator
# - Visualize the workflow (tasks and their connections)
# - Build and run a workflow with input parameters
# - Inspect results and the full provenance

# %%
# What’s Next
# -----------
#
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Concepts <../concept/index.rst>`__                     | A brief introduction of ``WorkGraph``’s main concepts. |
# |                                                         |                                                        |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `HowTo <../howto/index.rst>`__                          | Advanced topics and tips, e.g., flow control using     |
# |                                                         | ``if``, ``for``, ``while``, and ``context``.           |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Tutorials <../tutorial/index.rst>`__                   | Real-world examples in computational materials         |
# |                                                         | science and more.                                      |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Migration <migration_from_aiida_core/index.rst>`__     | Migration guide for users interested in transferring   |
# |                                                         | aiida-core processes (e.g. ``WorkChain``, ``CalcJob``) |
# |                                                         | to ``WorkGraph``.                                      |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
#
# .. _AiiDA documentation: https://aiida.readthedocs.io/en/stable/
