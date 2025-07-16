"""
=========================
Control flow in WorkGraph
=========================
"""
# %%
# Introduction
# ============
# In this how-to we show you how you can implement the ``while`` loop and ``if`` conditional flow control elements in ``WorkGraph``.
# You can implement both using the context manager approach, as well as using the ``graph_builder`` for the ``if`` conditional (``while`` is currently not supported by the ``graph_builder``).
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
# While
# =====

#%% 
# Workflow description
# --------------------
# Suppose we want to implement a workflow with the following tasks:


@task
def compare(x, y):
    return x < y


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# And the following workflow logic expressed in Python:
#
# .. code-block:: python
#
#    n = add(1, 1)
#    while compare(n, 8):
#        n = add(n, 1)
#        n = multiply(n, 2)
#    result = add(n, 1)

# %%
# Context manager
# ---------------
# Unlike regular tasks, the *while zone* lacks data input and output sockets.
# Tasks outside the zone can directly link to those inside, facilitating a workflow integration.
# We have the option to specify the maximum number of iterations to prevent an infinite loop.

with WorkGraph("while_context_manager") as wg:
    # set a context variable before running.
    outputs1 = add(x=1, y=1)
    wg.ctx.n = outputs1.result
    # We need to use context or static variables here
    outputs2 = compare(x=wg.ctx.n, y=8)
    # Because context variables do not follow the dataflow programming paradigm,
    # we have to state the dependency of the context variable explicitely.
    # The syntax below reads "`should_run` associated task waits for its executaion till `outputs1` has been successfully set its value"
    outputs2 << outputs1
    with While(outputs2.result, max_iterations=10):
        outputs3 = add(x=wg.ctx.n, y=1)
        outputs4 = multiply(x=outputs3.result, y=2)
        # We need to update the context variable since it is accessed in should_run task
        wg.ctx.n = outputs4.result
    wg.outputs.result = add(x=outputs4.result, y=1).result

wg.run()
print("State of WorkGraph:   {}".format(wg.state))
print("Result            :   {}".format(wg.outputs.result.value))
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view the ``While`` context manager is depicted as a *zone*, containing all its child tasks.
# This zone simplifies the visualization of the loop structure as it separates the logic executed within the loop from the one outside.
# Notice that the cyclic links around the context variable ``n`` since it is reused each iteration.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~

generate_node_graph(wg.pk)

# %%
# If
# ==

# %%
# Workflow description
# --------------------
#
# Now, using the same ``add`` and ``multiply`` tasks as defined for the ``while`` workflow, suppose you have the following workflow logic

# %%
# And the following workflow logic expressed in Python:
#
# .. code-block:: python
#
#    # First addition
#    result = add(1, 1)
#
#    # Conditionally execute addition or multiplication
#    if result < 0:
#        result = add(result, 2)
#    else:
#        result = multiply(result, 3)
#
#    # Last addition
#    result = add(result, 1)

# %%
# Context manager
# ---------------
#
# WorkGraph provides the ``If`` context manager that allows you to define conditional logic in your workflow.
# The ``If`` block encapsulates all its child tasks, which are executed based on the defined conditions.

with WorkGraph("if_context_manager") as wg:
    result = add(x=1, y=1).result
    with If(result < 0):
        r1 = add(x=result, y=2).result
        wg.ctx.result = r1
    with If(result >= 0):
        r2 = multiply(x=result, y=3).result
        wg.ctx.result = r2

    # -----------------------------------------
    result = add(x=wg.ctx.result, y=1).result
    r1 >> result
    r2 >> result
    wg.outputs.result = result

wg.run()
print("State of WorkGraph:   {}".format(wg.state))
print("Result            :   {}".format(wg.outputs.result.value))
# import ipdb; ipdb.set_trace()
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

generate_node_graph(wg.pk)

# %%
# Submit the WorkGraph
# ~~~~~~~~~~~~~~~~~~~~

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result:             {result.value}")

# %%
# Finally, after the WG has finished, we generate the node (provenance) graph from the AiiDA process, where we can see
# that the result of the ``op_lt`` (larger than) comparison is ``False``, while for the ``og_ge`` (greater or equal)
# comparison it is ``True``, meaning that the branch with the intermediate multiplication was executed.

generate_node_graph(wg.pk)


# %%
# Graph builder
# -------------
#
# The ``graph_builder`` decorator is used for creating a dynamic ``WorkGraph`` during runtime based on input values (see `this
# section <https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/graph_builder.html>`_).
#
# This method differs significantly from the ``If`` context manager:
#
# - **Visibility**: In the GUI, only the ``graph_builder`` task is visible before execution, while for the ``If``,
#   both branches were shown
# - **Dynamic Generation**: Upon running, it generates the WorkGraph dynmically, allowing for complex conditional logic and flow adjustments based on runtime data.

# Create a WorkGraph which is dynamically generated based on the input
# then we output the result as a graph-level output
@task.graph_builder(outputs=["result"])
def add_multiply_if(x, y):
    with WorkGraph("add_multiply_graph_builder") as wg:
        if x.value < 0:
            add1 = add(x=x, y=y)
            # export the result of add1 to the graph-level outputs
            wg.outputs.result = add1.outputs.result
        else:
            multiply1 = multiply(x=x, y=y)
            # export the result of multiply1 to the graph-level outputs
            wg.outputs.result = multiply1.outputs.result
    return wg

# 
with WorkGraph("if_graph_builer") as wg:
    add1 = wg.add_task(add, name="add1", x=1, y=1)
    add_multiply_if1 = wg.add_task(
        add_multiply_if, name="add_multiply_if1", x=add1.outputs.result, y=3
    )
    add1 = wg.add_task(add, name="add2", x=add_multiply_if1.outputs.result, y=1)

# We export the workgraph to an html file so that it can be visualized in a browser
wg.to_html()

# Comment out the following line to visualize the workgraph in jupyter-notebook
# wg

# %%
# In the GUI, we only see the ``add_multiply_if1`` task. When this task run, it will generate a ``WorkGraph`` based on the input value. This is different from the ``If`` task, in which we see all tasks before the WorkGraph run.


# %%
# Submit the WorkGraph and check the results
# ------------------------------------------

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.tasks.add2.outputs.result.value}")

# %%
# Generate node graph from the AiiDA process,and we can see that the ``multiply`` task is executed.

generate_node_graph(wg.pk)

# %%
# Using the If Task
# =================
#
# Internally, the ``If`` context manager is implemented using the ``If`` ``Task`` from the WorkGraph library.
# In the WorkGraph user interface, the ``If`` ``Task`` is visually represented as an "If Zone".
# This zone encapsulates all its child tasks, which are executed based on the defined conditions.
#
# - **conditions**: The If Zone includes a ``conditions`` socket, which determines when the tasks inside the zone should be executed.
# - **invert_condition**: If this input is True, it will invert the conditions.
# - **Task Linking**: Tasks located outside the If Zone can be directly linked to tasks within the zone, allowing for dynamic workflow adjustments based on conditional outcomes.
#
# Here is an example of how to add an ``If`` ``Task`` to a WorkGraph:
#
# .. code-block:: python
#
#     if_task = wg.add_task(
#        "workgraph.if_zone",
#        name="if_false",
#        conditions=condition1.outputs["result"],
#        invert_condition=True
#     )

# %%
# Adding tasks to the If Zone
# ---------------------------
# We can add tasks to the ``If`` zone using the ``children`` attribute.
#
# .. code-block:: python
#
#     # add task1 and task2 to the if zone
#     if_task.children.add(["task1", "task2"])


# %%
# Summary
# =======
#
# The ``If`` provides a visual and structured approach to managing conditional tasks within a defined zone. In contrast,
# the ``graph_builder`` decorator offers flexibility by dynamically generating the workflow based on runtime inputs,
# suitable for complex and adaptive process flows.

# %%
# Further reading
# ---------------
# Similarly, other the control logic can implemented, see `Use if condition  <../if.ipynb>`_ for if logic and `How to run tasks in parallel <parallel.py>`_ for map logic.
