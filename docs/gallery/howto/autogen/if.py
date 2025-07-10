"""
=====================================
Flow control: Using ``if`` conditions
=====================================
"""

# %%
# Introduction
# ============
#
# This tutorial provides a step-by-step guide on how to implement conditional logic in WorkGraph using three different
# methods:
#
# 1. **If context manager**
# 2. **graph_builder Decorator**
# 3. **Programmatic approach using the If Task**
#
# For simple cases, we recommend option 1), the ``If`` context manager approach, while option 2), using the Graph
# Builder provides additional advantages, see :doc:`graph_builder`.
# Finally, option 3) uses the ``If`` ``Task`` directly without the context manager. This approach requires a lot of
# boilerplate code and is generally not recommended. However, it presents the programmatic method if you want to
# construct the ``If`` flow control element using a node-graph programming approach.
#

# %%
# Using the ``If`` context manager
# ================================
#
# WorkGraph provides the ``If`` context manager that allows you to define conditional logic in your workflow. The ``If``
# block encapsulates all its child tasks, which are executed based on the defined conditions.

# %%
# Example
# -------
#
# Suppose you have the following Python workflow:


def add(x, y):
    return x + y


def multiply(x, y):
    return x * y


# step 1
result = add(1, 1)  # First addition
# step 2
if result < 0:
    result = add(result, 2)  # Conditionally execute addition
else:
    result = multiply(result, 2)  # Or multiplication
# step 3
result = add(result, 1)  # Last addition

print("Result is", result)

# %%
# To convert this workflow into a WorkGraph, we first convert the Python functions, ``add`` and ``multiply`` to WG tasks, as
# usual:

from aiida_workgraph import task


@task.calcfunction
def add(x, y):
    return x + y


@task.calcfunction
def multiply(x, y):
    return x * y


# %%
# To then define define the conditional logic for the workflow, we use the ``If`` context manager provided by
# AiiDA WorkGraph.
# Note that to use the ``If`` context manager, we also need to construct our top-level WorkGraph via the context manager
# approach. More information on this is provided in :doc:`../concept/workgraph`.

from aiida_workgraph import WorkGraph, If

with WorkGraph("if_context") as wg:
    result = add(x=1, y=1)
    with If(result < 0):
        wg.ctx.result = add(x=result, y=2)
    with If(result >= 0):
        wg.ctx.result = multiply(x=result, y=2)
    # -----------------------------------------
    result = add(x=wg.ctx.result, y=1)

# We export the workgraph to an html file so that it can be visualized in a browser
wg.to_html()

# Comment out the following line to visualize the workgraph in jupyter-notebook
# wg

# %%
# In the GUI, two boxes, or **Zone**\ s are being shown, one for each branch as defined by the ``If``, where each **If Zone**
# has a ``conditions`` socket, and their ``result``\ s are fed into the ``select`` **Zone**.

# %%
# Submit the WorkGraph and check the results
# ===========================================

from aiida import load_profile

load_profile()

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result:             {result.value}")

# %%
# Finally, after the WG has finished, we generate the node (provenance) graph from the AiiDA process, where we can see
# that the result of the ``op_lt`` (larger than) comparison is ``False``, while for the ``og_ge`` (greater or equal)
# comparison it is ``True``, meaning that the branch with the intermediate multiplication was executed.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Using the graph_builder Decorator
# ==================================
#
# The ``graph_builder`` decorator is used for creating a dynamic ``WorkGraph`` during runtime based on input values (see `this
# section <https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/graph_builder.html>`_).
#
# This method differs significantly from the ``If`` context manager:
#
# - **Visibility**: In the GUI, only the ``graph_builder`` task is visible before execution, while for the ``If``,
#   both branches were shown
# - **Dynamic Generation**: Upon running, it generates the WorkGraph dynamically, allowing for complex conditional logic and flow adjustments based on runtime data.

# Create a WorkGraph which is dynamically generated based on the input
# then we output the result as a graph-level output
# @task.graph_builder(outputs=[{"name": "result"}])
# def add_multiply_if(x, y):
#     wg = WorkGraph()
#     if x.value > 0:
#         add1 = wg.add_task(add, name="add1", x=x, y=y)
#         # export the result of add1 to the graph-lvel outputs
#         wg.outputs.result = add1.outputs.result
#     else:
#         multiply1 = wg.add_task(multiply, name="multiply1", x=x, y=y)
#         # export the result of multiply1 to the graph-lvel outputs
#         wg.outputs.result = multiply1.outputs.result
#     return wg

@task.graph_builder(outputs=[{"name": "result", "from": "ctx.result"}])
def add_multiply_if(x, y):
    wg = WorkGraph()
    if x.value > 0:
        add1 = wg.add_task(add, name="add1", x=x, y=y)
        # Use the correct socket name 'result', not 'sum'
        wg.update_ctx({"result": add1.outputs.result})
    else:
        multiply1 = wg.add_task(multiply, name="multiply1", x=x, y=y)
        # Use the correct socket name 'result', not 'sum'
        wg.update_ctx({"result": multiply1.outputs.result})
    return wg

# %%
# Create the workflow
# ===================

from aiida_workgraph import WorkGraph

wg = WorkGraph("if_graph_builer")
add1 = wg.add_task(add, name="add1", x=1, y=1)
add_multiply_if1 = wg.add_task(
    add_multiply_if, name="add_multiply_if1", x=add1.outputs.result, y=2
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
# ===========================================

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.tasks.add2.outputs.result.value}")
import ipdb; ipdb.set_trace()

# %%
# Generate node graph from the AiiDA process,and we can see that the ``multiply`` task is executed.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# If Task
# =======
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
# ===========================
# We can add tasks to the ``If`` zone using the ``children`` attribute.
#
# .. code-block:: python
#
#     # add task1 and task2 to the if zone
#     if_task.children.add(["task1", "task2"])

# %%
# Creating the workflow
# =====================
# To construct the workflow, we'll utilize the built-in ``If`` and ``Select`` tasks from the Workgraph library. The ``Select`` task enables us to choose between two data sources based on a specified condition.
#
# The ``Select`` task has the following inputs:
#
#    - **condition**: Provide the condition that dictates the selection between ``true`` and ``false`` outputs.
#    - **true**: Specify the output to be used if the condition evaluates to ``true``.
#    - **false**: Define the output for when the condition evaluates to ``false``.

from aiida_workgraph import WorkGraph

with WorkGraph("if_task") as wg:
    condition = add(x=1, y=1)

    if_true_zone = wg.add_task(
        "workgraph.if_zone", name="if_true", conditions=condition
    )
    add2 = if_true_zone.add_task(add, name="add2", x=condition, y=2)

    if_false_zone = wg.add_task(
        "workgraph.if_zone",
        name="if_false",
        conditions=condition,
        invert_condition=True,
    )
    multiply1 = if_false_zone.add_task(multiply, name="multiply1", x=condition, y=2)
    # ---------------------------------------------------------------------
    select1 = wg.add_task(
        "workgraph.select",
        name="select1",
        true=add2.outputs["result"],
        false=multiply1.outputs["result"],
        condition=condition,
    )
    add3 = wg.add_task(add, name="add3", x=select1.outputs["result"], y=1)

# We export the workgraph to an html file so that it can be visualized in a browser
wg.to_html()

# Comment out the following line to visualize the workgraph in jupyter-notebook
# wg

# %%
# Summary
# =======
#
# The ``If`` provides a visual and structured approach to managing conditional tasks within a defined zone. In contrast,
# the ``graph_builder`` decorator offers flexibility by dynamically generating the workflow based on runtime inputs,
# suitable for complex and adaptive process flows.
