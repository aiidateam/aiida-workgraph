"""
Control task execution order
============================
.. _task_execution_order:
"""

# %%
# Introduction
# ------------
#
# In `aiida-workgraph`, a task typically runs as soon as its inputs are available.
# However, you sometimes need to control the execution order more precisely, forcing
# a task to wait for another unrelated task to complete. This is useful when a task
# relies on data that is not passed directly as an input but is, for instance,
# updated in a shared context.
#
# .. note::
#
#    For dependencies between different WorkGraphs, use the ``monitor`` task.
#    Please refer to the :doc:`monitor <./monitor>` HowTo for more details.

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel('REPORT')
# sphinx_gallery_end_ignore

from aiida import load_profile
from aiida_workgraph import group, task

load_profile()


# %%
# Dependency operators ``>>`` and ``<<``
# --------------------------------------
#
# The simplest way to create an explicit dependency is with the ``>>`` and ``<<``
# operators. They let you enforce that one task must finish before another begins,
# even if there's no data link between them.
#
# * ``task_A >> task_B`` means **"task A, THEN task B"**.
# * ``task_B << task_A`` means **"task B WAITS ON task A"** - equivalent to ``task_A >> task_B``.
#
# .. note::
#
#    These operators can be used to create a chain of dependencies.
#
#    .. code-block:: python
#
#       # task_A THEN task_B THEN task_C
#       task_A >> task_B >> task_C
#
#    For many dependencies, the ``group`` utility can be used:
#
#    .. code-block:: python
#
#       from aiida_workgraph import group
#
#       # tasks A AND B (in parallel), THEN task C
#       group(task_A, task_B) >> task_C


# %%
# Simple example
# --------------
#
# Let's define a workflow where an ``add`` task must wait for two ``multiply``
# tasks to finish first.


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# Next, we define a WorkGraph in which a group of ``multiply`` tasks run in parallel,
# then (``>>``) an ``add`` task executes, returning its result.


@task.graph
def wait_graph(x, y):
    return (
        group(
            multiply(x=y, y=2),
            multiply(x=y, y=3),
        )
        >> add(x=x, y=1).result  # `add` waits for both `multiply` tasks to complete
    )


wg = wait_graph.build(x=1, y=2)
wg.run()

# %%
# By checking the ``REPORT`` logs from AiiDA, you will see that both ``multiply``
# tasks complete before the ``add`` task begins, as specified.


# %%
# Conclusion
# ----------
#
# You now know two effective methods to control task execution order in ``aiida-workgraph``:
#
# * Use the wait operators (``>>``, ``<<``) for simple, direct dependencies
#   between individual tasks.
# * Use a ``Zone`` to group tasks that should be treated as a single execution
#   block, simplifying complex dependency chains.
#
# These tools give you fine-grained control over your workflows, ensuring your
# calculations run in the correct sequence.

# sphinx_gallery_start_ignore
set_aiida_loglevel('ERROR')
# sphinx_gallery_end_ignore
