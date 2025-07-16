"""
=========================
Control flow in WorkGraph
=========================
"""
# %%
# Introduction
# ============
# In this how-to we show you how you can implement the ``while`` loop and ``if`` conditional flow control elements in ``WorkGraph``.
# You can implement both using the context manager approach, as well as using the ``@task.graph`` decorator for the ``if`` conditional (``while`` is currently not supported).
# So let's dive right into it!

# %%
# Setting up the AiiDA environment
# --------------------------------
# First, we need to set up our AiiDA environment, importing the necessary entities from ``aiida-workgraph`` and loading the AiiDA profile.

from aiida_workgraph import WorkGraph, task, If, While
from aiida_workgraph.utils import generate_node_graph
from aiida import load_profile

load_profile()


# %%
# If
# ==

# %%
# Workflow description
# --------------------
#
# Now, let's see how we can convert the following arithmetic workflow containing an ``if`` conditional into a WorkGraph:
#
# .. code-block:: python
#
#    # First addition
#    result = 1 + 1
#
#    # Conditionally execute addition or multiplication
#    if result < 0:
#        result = result + 2
#    else:
#        result = result * 3
#
#    # Last addition
#    result = result + 1

# %%
# We now define the relevant arithmetic operations as WorkGraph tasks.
# They will present the processes executed in the workflow and their provenance will be tracked.

@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y

# %%
# Context manager
# ---------------
#
# For this, WorkGraph also provides an ``If`` context manager, which allows you to define conditional logic in your workflow.
# The ``If`` block encapsulates all its child tasks, which are executed based on the defined conditions.

with WorkGraph("if_task") as wg:
    outputs1 = add(x=1, y=1)
    with If(outputs1.result < 0):
        outputs2 = add(x=outputs1.result, y=2)
        wg.ctx.result = outputs2.result
    with If(outputs1.result >= 0):
        outputs3 = multiply(x=outputs1.result, y=3)
        wg.ctx.result = outputs3.result

    outputs4 = add(x=wg.ctx.result, y=1)
    # Again, we have to use the ``>>`` syntax here to explicitly tell WG to wait with the execution of the last task,
    # until the previous two tasks have finished.
    outputs2 >> outputs4
    outputs3 >> outputs4
    # In principle, as we are dealing with python functions that are run in a blocking manner here, the example would
    # also work without ``>>``. However, if the tasks would be submitted to the deamon in a non-blocking fashion, the
    # explicit waiting enforced by ``>>`` is required.

    # Finally, we set the result of the last task as the global workflow output result
    wg.outputs.result = outputs4.result

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.outputs.result.value}")

assert wg.outputs.result.value == 7

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view, two boxes, or ``zone``s are being shown, one for each branch as defined by the ``If``
# Here, each ``If`` ``zone`` has a ``conditions`` socket, and the ``result``s are fed into the ``Select`` ``zone``.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# Finally, after the WG has finished, we generate the node (provenance) graph from the AiiDA process, where we can see
# that the result of the ``op_lt`` (larger than) comparison is ``False``, while for the ``og_ge`` (greater or equal)
# comparison it is ``True``, meaning that the branch with the intermediate multiplication was executed.

generate_node_graph(wg.pk)

# %%
# Graph decorator
# ---------------
#
# The ``graph`` decorator is used for creating a dynamic ``WorkGraph`` during runtime based on input values (see `this
# section <https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/graph.html>`_).
#
# This method differs significantly from the ``If`` context manager:
#
# - **Visibility**: In the GUI, only the ``graph`` task is visible before execution, while for the ``If``,
#   both branches were shown
# - **Dynamic Generation**: Upon running, it generates the WorkGraph dynamically, allowing for complex conditional logic and flow adjustments based on runtime data.

# Create a WorkGraph which is dynamically generated based on the input
# then we output the result as a graph-level output
@task.graph()
def add_multiply_if(x, y, z):
    if x.value < 0:
        return add(x=x, y=y).result
    else:
        return multiply(x=x, y=z).result

with WorkGraph("if_graph") as wg:
    add1 = wg.add_task(add, name="add1", x=1, y=1)
    add_multiply_if1 = wg.add_task(
        add_multiply_if, name="add_multiply_if1", x=add1.outputs.result, y=2, z=3
    )
    add1 = wg.add_task(add, name="add2", x=add_multiply_if1.outputs.result, y=1)

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.tasks.add2.outputs.result.value}")

# %%
# Workflow view
# ~~~~~~~~~~~~~
#
# TODO: Continue here
# In the GUI, we only see the ``add_multiply_if1`` task. When this task run, it will generate a ``WorkGraph`` based on the input value. This is different from the ``If`` task, in which we see all tasks before the WorkGraph run.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# Generate node graph from the AiiDA process,and we can see that the ``multiply`` task is executed.

generate_node_graph(wg.pk)

# %%
# While
# =====

# %%
# Workflow description
# --------------------
# Suppose we want to implement a workflow with the following logic:

# .. code-block:: python
#
#    n = 1 + 1
#    while n < 8:
#        n = n + 1
#        n = n * 2
#        print(n)
#    result = n + 1

# %%
# We have to define our comparison operation as a task, as well
# PRCOMMENT: Why??

@task
def compare_lt(x, y):
    return x < y

# %%
# Context manager
# ---------------
# Now, we can construct the WorkGraph using the ``While`` context manager,
# which effectively creates a ``while zone`` (as seen in the workflow view below).
# Unlike regular tasks, the *while zone* lacks data input and output sockets.
# Tasks outside the zone can directly link to those inside, facilitating workflow integration.
# We have the option to specify the maximum number of iterations to prevent an infinite loop.

with WorkGraph("while_context_manager") as wg:
    initial_n = add(x=1, y=1).result
    wg.ctx.n = initial_n
    should_run = compare_lt(x=wg.ctx.n, y=8).result
    should_run << initial_n
    with While(should_run, max_iterations=10):
        n = add(x=wg.ctx.n, y=1).result
        n = multiply(x=n, y=2).result
        wg.ctx.n = n

    wg.outputs.result = add(x=n, y=1).result

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.outputs.result.value}")
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15


# %% Notes on the example before

# set a context variable before running.
# We need to use context or static variables here
# Because context variables do not follow the dataflow programming paradigm,
# we have to state the dependency of the context variable explicitly.
# The syntax below reads "``should_run`` associated task waits for its
# execution till ``outputs1`` has been successfully set its value"

# We need to update the context variable since it is accessed in should_run task

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view the ``While`` context manager is depicted as a *while zone*, containing all its child tasks.
# This zone simplifies the visualization of the loop structure as it separates the logic executed within the loop from the one outside.
# Notice that the cyclic links around the context variable ``n`` since it is reused each iteration.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~

generate_node_graph(wg.pk)