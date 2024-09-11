"""
==========================================================
Use Graph Builder to create a nested and dynamic workflows
==========================================================
"""

# %%
# Nested workflows
# ================
# The `Graph Builder` allow user to create nested workflows from an input.

# Load the AiiDA profile.
from aiida import load_profile

load_profile()


# %%
# Example
# -------
# Suppose we want a WorkGraph which includes another WorkGraph`(x+y)*z` inside it.
# We can actually add a WorkGraph to another WorkGraph

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


def add_multiply(x, y, z):
    wg = WorkGraph()
    wg.add_task(add, name="add", x=x, y=y)
    wg.add_task(multiply, name="multiply", x=z)
    wg.add_link(wg.tasks["add"].outputs[0], wg.tasks["multiply"].inputs["y"])
    return wg


wg = WorkGraph("nested_workgraph")
add_multiply1 = wg.add_task(add_multiply(x=Int(2), y=Int(3), z=Int(4)))
wg.to_html()

# %%
# Run the workgraph

wg.run()

# %%
# However linking the two WorkGraphs will not work

wg = WorkGraph("nested_workgraph")
add_multiply1 = wg.add_task(add_multiply(x=Int(2), y=Int(3), z=Int(4)))

try:
    wg.add_task(
        add_multiply(x=add_multiply1.outputs["multiply.result"], y=Int(3), z=Int(4))
    )
except Exception as err:
    print(err)

# %%
# For that use case we need to use the graph builder

# Create a graph builder function
# -------------------------------
# We add `task.graph_builder` decorator to a function to define a graph builder
# function. The function constructs a WorkGraph based on the inputs, and returns
# it at the end.
#
#
# Expose outputs
# --------------
# We can expose the outputs of the tasks as the outputs of the WorkGraph:
#
# .. code:: python
#
#     @task.graph_builder(outputs = [{"name": "multiply", "from": "multiply.result"}])
#
# This will expose the `result` output of the `multiply` task as the `multiply` output of the WorkGraph.
#
#


# We use task.graph_builder decorator and expose the output of the "multiply"
# task as the output of the graph builder function.
@task.graph_builder(outputs=[{"name": "multiply", "from": "multiply.result"}])
def add_multiply(x, y, z):
    # Create a WorkGraph
    wg = WorkGraph()
    wg.add_task(add, name="add", x=x, y=y)
    wg.add_task(multiply, name="multiply", x=z)
    wg.add_link(wg.tasks["add"].outputs[0], wg.tasks["multiply"].inputs["y"])
    # Don't forget to return the `wg`
    return wg


# %%
# Create nested workflow
# ----------------------
# We can use the graph builder task inside another WorkGraph to create nested workflow. Here is an example:


wg = WorkGraph("test_graph_build")
# create a task using the graph builder
add_multiply1 = wg.add_task(add_multiply, x=Int(2), y=Int(3), z=Int(4))
add_multiply2 = wg.add_task(add_multiply, x=Int(2), y=Int(3))
# link the output of a task to the input of another task
wg.add_link(add_multiply1.outputs[0], add_multiply2.inputs["z"])
wg.submit(wait=True)
assert add_multiply2.outputs[0].value == 100
wg.to_html()


# %%
# Generate node graph from the AiiDA process,and we can see that the `multiply` task is executed.


from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Looking at the process list we can also that multiple WorkGraphs have been submitted
#
# .. code-block:: bash
#
#     verdi process list -a


# %%
# Use the graph builder directly
# -------------------------------
# Of course, one can use the graph builder directly to create a WorkGraph. Here is an example:


wg = add_multiply(2, 3, 4)
wg.submit(wait=True)


# %%
# Create a Task from the workgraph (Experimental)
# -----------------------------------------------
# One can create a Task from a WorkGraph directly.


from aiida_workgraph import WorkGraph

wg1 = WorkGraph()
# Note, one can not set the inputs values here using AiiDA data types
wg1.add_task(add, name="add")
wg1.add_task(multiply, name="multiply")
wg1.add_link(wg1.tasks["add"].outputs[0], wg1.tasks["multiply"].inputs["y"])


# %%
# Then we can use this `wg1` inside a WorkGraph. The inputs and outputs of all tasks in `wg1` will be exposed as the inputs and outputs of the task.


wg2 = WorkGraph("test_graph_build")
# create a task using the graph builder
add_multiply1 = wg2.add_task(wg1, name="add_multiply1")
add_multiply2 = wg2.add_task(wg1, name="add_multiply2")
# link the output of a task to the input of another task
wg2.add_link(add_multiply1.outputs["multiply.result"], add_multiply2.inputs["add.x"])


# Create a task using the WorkGraph
print("Inputs:")
for input in add_multiply1.inputs:
    print(f"  - {input.name}")
print("Outputs:")
for output in add_multiply1.outputs:
    print(f"  - {output.name}")

wg2.to_html()


# %%
# Prepare the inputs and run the WorkGraph:


from aiida.orm import Int

add_multiply1.set({"add": {"x": Int(2), "y": Int(3)}, "multiply": {"x": Int(4)}})
add_multiply2.set({"add": {"y": Int(4)}, "multiply": {"x": Int(3)}})

wg2.submit(wait=True)


# %%
# Generate node graph from the AiiDA process

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg2.pk)


# %%
# More usage (like `if` and `while`) of graph builder will be shown in the following tutorials.

# %%
# Dynamic workflows
# =================
# The `Graph Builder` also allows us to create dynamic workflows that can change depending on the input.


# %%
# Example for loop
# ----------------
# In this example we will create a dynamic number of tasks as specified in the
# input of the WorkGraph.


@task.calcfunction()
def add_one(x):
    return x + 1


@task.graph_builder(outputs=[{"name": "result", "from": "context.task_out"}])
def for_loop(nb_iterations: Int):
    wg = WorkGraph()
    for i in range(nb_iterations):
        task = wg.add_task(add_one, x=i)

    # We cannot refer to a specific task as output in the graph builder decorator
    # as in the examples below as the name of the last task depends on the input.
    # Therefore we use the context to not directly refer to the name but the last
    # task object that was created. The context can then be reffered in the outputs
    # of the graph builder decorator.
    task.set_context(
        {"result": "task_out"}
    )  # put result of the task to the context under the name task_out
    # If want to know more about the usage of the context please refer to the
    # context howto in the documentation
    return wg


wg = WorkGraph()
task = wg.add_task(for_loop, nb_iterations=Int(2))
wg.to_html()

# %%
# Running the workgraph.

wg.submit(wait=True)
print("Output of last task", task.outputs["result"].value)  # 1 + 1 result

# %%
# Plotting provenance

generate_node_graph(wg.pk)


# %%
# Example if-then-else
# --------------------
# Suppose we want to run a different task depending on the input. We run the
# add_one task if the number is below 2 otherwise we run a modulo 2
# task.

from aiida_workgraph import task, WorkGraph
from aiida.orm import Int


@task.calcfunction()
def modulo_two(x):
    return x % 2


@task.graph_builder(outputs=[{"name": "result", "from": "context.task_out"}])
def if_then_else(i: Int):
    wg = WorkGraph()
    if i.value < 2:
        task = wg.add_task(add_one, x=i)
    else:
        task = wg.add_task(modulo_two, x=i)

    # same concept as before, please read the for loop example for explanation
    task.set_context({"result": "task_out"})
    return wg


wg = WorkGraph()
task1 = wg.add_task(if_then_else, i=Int(1))
task2 = wg.add_task(if_then_else, i=task1.outputs["result"])
wg.to_html()

# %%
# Running the workgraph.

wg.submit(wait=True)
print("Output of first task", task1.outputs["result"].value)  # 1 + 1 result
print("Output of second task", task2.outputs["result"].value)  # 2 % 2 result

# %%
# Plotting provenance

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)
