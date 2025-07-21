"""
=================================
Control task execution order
=================================

In ``aiida-workgraph``, a task typically runs as soon as its inputs are available.
However, you sometimes need to control the execution order more precisely, forcing
a task to wait for another unrelated task to complete. This is useful when a task
relies on data that is not passed directly as an input but is, for instance,
updated in a shared context.

This tutorial demonstrates two powerful features for managing these dependencies:

* The wait operators: ``>>`` and ``<<``
* Task grouping with ``Zone``

.. note::

   For dependencies between different WorkGraphs, use the ``monitor`` task.
   Please refer to the :doc:`monitor <./monitor>` tutorial for more details.

"""
# %%
# Initial setup
# -------------
# First, let's set up the environment and adjust the AiiDA log level to ``REPORT`` so we can observe the execution order of the tasks.

from aiida_workgraph.utils.logging import set_aiida_loglevel
from aiida_workgraph import task, WorkGraph, Zone
from aiida import load_profile

set_aiida_loglevel("REPORT")

load_profile()


# %%
# Explicit dependencies with wait operators (``>>`` and ``<<``)
# ==========================================================
#
# The simplest way to create an explicit dependency is with the ``>>`` and ``<<``
# operators. They let you enforce that one task must finish before another begins,
# even if there's no data link between them.
#
# * ``task_A >> task_B`` means **"task B waits for task A"**.
# * ``task_B << task_A`` means the exact same thing.
#
# Example
# -------
# Let's define a workflow where an ``add`` task must wait for two ``multiply``
# tasks to finish first.


@task()
def add(x, y):
    return x + y


@task()
def multiply(x, y):
    return x * y


@task.graph()
def wait_graph(x, y):
    add1_outputs = add(x=x, y=1)
    multiply1_outputs = multiply(x=y, y=2)
    multiply2_outputs = multiply(x=y, y=3)

    # Make the 'add_task' wait for both 'multiply' tasks to complete
    multiply1_outputs >> add1_outputs
    multiply2_outputs >> add1_outputs

    return add1_outputs.result


# %%
# .. note::
#
#    For many dependencies, the ``group`` utility can be used:
#
#    .. code-block:: python
#
#       from aiida_workgraph.collection import group
#       group(multiply1_outputs, multiply2_outputs) >> add1_outputs

# Build and run the WorkGraph
wg = wait_graph.build_graph(x=1, y=2)
wg.run()

# %%
# By checking the ``REPORT`` logs from AiiDA, you will see that both ``multiply``
# tasks complete before the ``add`` task begins, just as we specified.


# %%
# Grouping Dependencies with ``Zone``
# ===================================
#
# For more complex scenarios, you can group a set of tasks into a **Zone**.
# A ``Zone`` acts as a single unit for dependency management, governed by two rules:
#
# 1. **Entry Condition**: A ``Zone`` (and all tasks within it) will only start
#    after *all* tasks with links pointing *into* the ``Zone`` are finished.
# 2. **Exit Condition**: Any task that needs an output from *any* task inside the
#    ``Zone`` must wait for the entire `Zone` to complete.
#
# Example
# -------
# Here, we'll create a workflow where ``task3`` depends on ``task1``, but because
# ``task1`` is inside a ``Zone``, ``task3`` must wait for the whole group
# (``task1`` and ``task2``) to finish.

# Create a WorkGraph
with WorkGraph("zone_example") as wg:
    task0_outputs = add(x=1, y=1)

    # This Zone will only start after task1 is finished,
    # because task3 depends on its result.
    with Zone() as zone1:
        task1_outputs = add(x=1, y=1)
        task2_outputs = add(x=1, y=task0_outputs.result)

    # Task 4 will wait for the entire Zone to finish,
    # even though it only needs the result from task2.
    task3_outputs = add(x=1, y=task1_outputs.result)

# Generate an HTML file to visualize the graph
wg.to_html()

# %%
# Run the WorkGraph
wg.run()

# %%
# If you examine the logs, you'll see that ``task3``
# only starts after both ``task1`` and ``task2`` are complete, demonstrating the
# "all or nothing" behavior of a ``Zone``.


# %%
# Conclusion
# ==========
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

# Reset the log level to avoid verbose output in other examples
set_aiida_loglevel("ERROR")
