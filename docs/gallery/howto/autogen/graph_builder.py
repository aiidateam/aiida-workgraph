"""
==============================================
Graph Builder for nested and dynamic workflows
==============================================
"""

# %%
# Introduction
# ============
# In this example we learn how to create nested workflows by creating a task
# out of a WorkGraph. Further, we will learn how to do the same with the Graph
# Builder, a decorator that allows us to move the creation of the WorkGraph to
# runtime, so we can create dynamic workflows that change depending on the inputs.
# This is of particular interest for integrating for-loops and if-then-else
# logic into your workflow.

# Load the AiiDA profile.
from aiida import load_profile

load_profile()

# %%
# Nested workflows with WorkGraph
# ===============================
# We will discuss how to use WorkGraph's
# nested workflows. Suppose we want to reuse the WorkGraph computing `(x+y)*z`
# to perform the operation
#
# .. code-block:: Python
#
#    wg_out = (x+y)*z
#    out = (x+y)*wg_out
#
# We can integrate a WorkGraph to another WorkGraph by creating a task out of it.

from aiida_workgraph import task, WorkGraph
from aiida.orm import Int

# define add task
@task.calcfunction()
def add(x, y):
    return x + y


# define multiply task
@task.calcfunction()
def multiply(x, y):
    return x * y


def add_multiply(x=None, y=None, z=None):
    wg = WorkGraph()
    wg.add_task(add, name="add", x=x, y=y)
    wg.add_task(multiply, name="multiply", x=z)
    wg.add_link(wg.tasks.add.outputs.result, wg.tasks.multiply.inputs.y)
    return wg


wg = WorkGraph("nested_workgraph")
# Creating a task from the WorkGraph
add_multiply1 = wg.add_task(
    add_multiply(x=Int(2), y=Int(3), z=Int(4)), name="add_multiply1"
)
add_multiply2 = wg.add_task(add_multiply(x=Int(2), y=Int(3)), name="add_multiply2")
# link the output of a task to the input of another task
wg.add_link(add_multiply1.outputs.multiply.result, add_multiply2.inputs.multiply.x)
wg.to_html()

# %%
# The created WorkGraphTask behaves similarly as a normal WorkGraph would (and indeed actually has the associated
# `WorkGraph` attached as an attribute).That means we can access elements of the sub-WorkGraph, for instance, its tasks,
# inputs, etc., via:

print(wg.tasks.add_multiply1.tasks)
print(wg.tasks.add_multiply1.tasks.add.inputs.x)
# or
print(wg.tasks["add_multiply1"].tasks)
# and
print(wg.tasks.add_multiply1.inputs)
print(wg.tasks.add_multiply1.outputs)


# %%
# Finally, we run the workgraph

wg.submit(wait=True)
# (2+3)*4 = 20
# (2+3)*20 = 100
assert add_multiply2.outputs.multiply.result.value == 100

# %%
# And to generate the node graph from the AiiDA process

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Graph builder
# =============
# A much more powerful tool to create nested WorkGraphs is the Graph Builder.
# It is a decorator that we can add to a function that returns a WorkGraph
# similar as `add_multiply` to have more control what we want to expose to the
# user and to create dynamic workflows.


# %%
# Expose outputs
# --------------
# We add `task.graph_builder` decorator to a function to define a graph builder
# function. The function constructs a WorkGraph based on the inputs, and returns
# it at the end.
# We can expose the outputs of the tasks as the outputs of the WorkGraph:
#
# .. code-block:: python
#
#     @task.graph_builder(outputs = [{"name": "multiply", "from": "multiply.result"}])
#
# This will expose the `result` output of the `multiply` task as the `multiply` output of the WorkGraph.
#


# We use task.graph_builder decorator and expose the output of the "multiply"
# task as the output of the graph builder function.
@task.graph_builder(outputs=[{"name": "multiply", "from": "multiply.result"}])
def add_multiply(x, y, z):
    # Create a WorkGraph
    wg = WorkGraph()
    wg.add_task(add, name="add", x=x, y=y)
    wg.add_task(multiply, name="multiply", x=z)
    wg.add_link(wg.tasks.add.outputs[0], wg.tasks.multiply.inputs.y)
    # Don't forget to return the `wg`
    return wg


# %%
# Create nested workflow
# ----------------------
# We can use the graph builder task inside another WorkGraph to create nested
# workflow similar as with a regular WorkGraph.


wg = WorkGraph("test_graph_build")
# create a task using the graph builder, note the difference as the inputs
# are specified as ports here
add_multiply1 = wg.add_task(add_multiply, x=Int(2), y=Int(3), z=Int(4))
add_multiply2 = wg.add_task(add_multiply, x=Int(2), y=Int(3))
# link the output of a task to the input of another task
wg.add_link(add_multiply1.outputs[0], add_multiply2.inputs.z)
wg.submit(wait=True)
assert add_multiply2.outputs[0].value == 100
wg.to_html()


# %%
# Generate node graph from the AiiDA process,and we can see that the `multiply` task was executed.

generate_node_graph(wg.pk)

# %%
# Looking at the process list we can also that multiple WorkGraphs have been submitted.
# Please run this now in the terminal:
#
# .. code-block:: bash
#
#     verdi process list -a


# %%
# Use the graph builder directly
# ------------------------------
# Of course, one can use the graph builder directly to create a WorkGraph. Here is an example:

wg = add_multiply(2, 3, 4)
wg.submit(wait=True)


# %%
# More usage (like `if` and `while`) of graph builder will be shown in the following how-tos.

# %%
# Dynamic workflows
# -----------------
# The `Graph Builder` also allows us to create dynamic workflows that can change depending on the input.


# %%
# Example for loop
# ^^^^^^^^^^^^^^^^
# In this example we will create a dynamic number of tasks as specified in the
# input of the WorkGraph.


@task.calcfunction()
def add_one(x):
    return x + 1


@task.graph_builder(outputs=[{"name": "result", "from": "ctx.task_out"}])
def for_loop(nb_iterations: Int):
    wg = WorkGraph()
    for i in range(nb_iterations.value):
        task = wg.add_task(add_one, x=i)

    # We cannot refer to a specific task as output in the graph builder decorator
    # as in the examples before since the name of the last task depends on the input.
    # Remember that each task is always assigned unique name automatically.
    # Therefore we use the context to not directly refer to the name but the last
    # task object that was created. The context can then be referred in the outputs
    # of the graph builder decorator.

    # Put result of the task to the context under the name task_out
    task.update_ctx({"task_out": "result"})
    # If want to know more about the usage of the context please refer to the
    # context howto in the documentation
    return wg


wg = WorkGraph("Nested workflow: For")
loop_task = wg.add_task(for_loop, nb_iterations=Int(2))
wg.to_html()

# %%
# Running the workgraph.

wg.submit(wait=True)
print("Output of last task", loop_task.outputs.result.value)  # 1 + 1 result

# %%
# Plotting provenance

generate_node_graph(wg.pk)


# %%
# Example if-then-else
# ^^^^^^^^^^^^^^^^^^^^
# Suppose we want to run a different task depending on the input. We run the
# add_one task if the number is below 2 otherwise we run a modulo 2
# task.

from aiida_workgraph import task, WorkGraph
from aiida.orm import Int


@task.calcfunction()
def modulo_two(x):
    return x % 2


@task.graph_builder(outputs=[{"name": "result", "from": "ctx.task_out"}])
def if_then_else(i: Int):
    wg = WorkGraph()
    if i.value < 2:
        task = wg.add_task(add_one, x=i)
    else:
        task = wg.add_task(modulo_two, x=i)

    # same concept as before, please read the for loop example for explanation
    task.update_ctx({"task_out": "result"})
    return wg


wg = WorkGraph("Nested workflow: If")
task1 = wg.add_task(if_then_else, i=Int(1))
task2 = wg.add_task(if_then_else, i=task1.outputs.result)
wg.to_html()

# %%
# Running the workgraph.

wg.submit(wait=True)
print("Output of first task", task1.outputs.result.value)  # 1 + 1 result
print("Output of second task", task2.outputs.result.value)  # 2 % 2 result

# %%
# Plotting provenance

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)
