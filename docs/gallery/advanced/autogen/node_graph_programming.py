"""
=========================================================
Write workflows using the node-graph programming paradigm
=========================================================
"""

# %%
# This guide introduces the **node-graph programming** in `aiida-workgraph`, which provides an alternative approach for constructing workflows.
# Unlike Pythonic approach used elsewhere in the documentation, node-graph programming is the low-level approach where you build the graph piece by piece.
# You manually add each `task` and connect them with `links`.
# This method offers maximum control but is more verbose and is generally reserved for advanced use cases, like programmatically generating a graph's structure.
#
# We'll explore how to build workflows by adding individual tasks and linking their inputs and outputs.
# This includes simple sequential workflows, as well as more complex structures involving control flow like `if` conditions and `while` loops.
#
# First, we set up our AiiDA environment:

from aiida_workgraph import WorkGraph, task, spec
from aiida import load_profile

load_profile()

# %%
# Creating a Simple Workflow
# ==========================
#
# Let's define a few simple Python functions that will serve as our tasks.


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
# Now, building a workflow via the node-graph programming approach involves three core steps:
#
# 1.  **Instantiate an empty WorkGraph**: This is the container for your workflow.
# 2.  **Add tasks**: Incorporate the defined Python tasks into the ``WorkGraph``.
# 3.  **Link tasks**: Connect the outputs of one task to the inputs of another, defining
#     the data flow.
#
# Here's an example demonstrating an "add then multiply" workflow:

# 1. Create an empty WorkGraph
wg = WorkGraph("add_multiply_workflow")

# 2. Add tasks to the workgraph
add_task = wg.add_task(add, name="add1")
multiply_task = wg.add_task(multiply, name="multiply1")

# 3. Link the output of 'add1' to the 'x' input of 'multiply1'
wg.add_link(add_task.outputs.result, multiply_task.inputs.x)

# Define the graph-level outputs
wg.outputs.result = multiply_task.outputs.result

# Run the workflow with specific input values
wg.run(
    inputs={"add1": {"x": 2, "y": 3}, "multiply1": {"y": 4}},
)

print(f"State of WorkGraph: {wg.state}")
print(f"Result: {wg.outputs.result.value}")

# Visualize the workgraph
wg.to_html()

# %%
# Conditional logic with the `If Zone` task
# ==========================================
# `If` logic in `aiida-workgraph` is represented by an **If Zone**, which visually encapsulates child tasks that execute based on specific conditions.
#
# Key features of the If Zone:
#
# * **conditions socket**: Determines when tasks within the zone are executed.
# * **invert_condition**: If ``True``, reverses the outcome of the ``conditions``.
# * **Task Linking**: Tasks outside the If Zone can directly link to tasks inside, allowing dynamic workflow adjustments based on conditional outcomes.
#
# `aiida-workgraph` provides the built-in `workgraph.if_zone` task to create an `If Zone` and the `workgraph.select` task to choose between different data sources based on a condition.
#
# Let's build a workflow where the result of an initial `add` operation dictates whether a subsequent `add` or `multiply` operation is performed.

wg = WorkGraph("if_task_example")
add_task = wg.add_task(add, x=1, y=1)
# If the condition is true
if_true_zone = wg.add_task(
    "workgraph.if_zone", name="if_true", conditions=add_task.outputs.result
)
add2 = if_true_zone.add_task(
    add, name="add2", x=add_task.outputs.result, y=2
)  # 2 + 2 = 4
# If the condition is false
if_false_zone = wg.add_task(
    "workgraph.if_zone",
    name="if_false",
    conditions=add_task.outputs.result,
    invert_condition=True,
)
multiply1 = if_false_zone.add_task(
    multiply, name="multiply1", x=add_task.outputs.result, y=2
)  # 2 * 2 = 4
# Select the result based on the initial condition
select1 = wg.add_task(
    "workgraph.select",
    name="select1",
    true=add2.outputs["result"],
    false=multiply1.outputs["result"],
    condition=add_task.outputs.result,
)
# Add 1 to the selected result
add3 = wg.add_task(add, name="add3", x=select1.outputs["result"], y=1)  # 4 + 1 = 5
# Graph-level output
wg.outputs.result = add3.outputs.result

# Run the workflow
wg.run()

print(f"State of WorkGraph: {wg.state}")
print(f"Result: {wg.outputs.result.value}")
assert wg.outputs.result.value == 5

# Visualize the workgraph
wg.to_html()

# %%
# Loops with `While Zone` task
# =============================
# The `While` loop in `aiida-workgraph` functions similarly to programming `while` loops,
# repeatedly executing a set of tasks as long as a specified condition remains `True`.
# This is handled by the `workgraph.while_zone` task.

wg = WorkGraph("while_task_example")

# Initialize 'n' with an initial value
initial_add_task = wg.add_task(add, x=1, y=1)  # n = 2
wg.ctx.n = initial_add_task.outputs.result

# Define the condition for the while loop: n < 8
# Here, we use the `compare` task as defined above
condition_task = wg.add_task(compare, x=wg.ctx.n, y=8)
# Ensure the condition task waits for the initial_add_task to complete
condition_task.waiting_on.add(initial_add_task)

# Start the While Zone
while_task = wg.add_task(
    "workgraph.while_zone", max_iterations=10, conditions=condition_task.outputs.result
)

# Tasks within the while loop
# First, add 1 to n
add_task_in_loop = while_task.add_task(add, x=wg.ctx.n, y=1)
# Then, multiply the result by 2
multiply_task_in_loop = while_task.add_task(
    multiply, x=add_task_in_loop.outputs.result, y=2
)
# Update 'n' for the next iteration of the loop
wg.ctx.n = multiply_task_in_loop.outputs.result

# After the loop, add 1 to the final 'n'
final_add_task = wg.add_task(add, x=multiply_task_in_loop.outputs.result, y=1)
wg.outputs.result = final_add_task.outputs.result

# Run the workflow
wg.run()

print(f"State of WorkGraph: {wg.state}")
print(f"Result: {wg.outputs.result.value}")

assert wg.outputs.result.value == 15

# Visualize the workgraphs
wg.to_html()

# %%
# Mapping operations with `Map Zone` task
# =========================================
# .. warning::
#   **This feature is experimental.** The API for ``Map`` zone is subject to change in future releases. We welcome your feedback on its functionality.
#
# The `Map` task in `aiida-workgraph` allows you to apply a function or a set of tasks to each item in a dictionary,
# similar to Python's built-in `map()` function.
# This is particularly useful for parallelizing operations over a dataset.
#
# First, let's define a task that generates a dictionary of data. Notice the `outputs` decorator,
# which indicates that `result` is a dynamic output and will be a namespace.
from typing import Any


@task
def generate_data(N) -> spec.namespace(result=spec.dynamic(Any)):
    """Generates a dictionary with N items."""
    data = {f"item_{i}": i for i in range(N)}
    return {"result": data}


@task
def calc_sum(**kwargs):
    """Calculates the sum of all keyword arguments' values."""
    return sum(kwargs.values())


# %%
# To use the ``Map`` task, you define a ``source`` which is a dictionary.
# The tasks inside the `map_zone` will be executed for each `item` in the source,
# where `item` represents the value of each key-value pair.
#

wg = WorkGraph("map_task_example")

# Generate a dictionary of data with 4 items (0, 1, 2, 3)
data_task = wg.add_task(generate_data, N=4)

# Create a Map Zone, with the source being the dictionary from generate_data
map_task = wg.add_task("workgraph.map_zone", source=data_task.outputs.result)

# Inside the Map Zone, add 1 to each item
add_task_in_map = map_task.add_task(add, x=map_task.item, y=1)

# After the Map Zone, sum all the results from the add_task_in_map
# The 'kwargs' input allows collecting all dynamic outputs from the mapped tasks.
sum_task = wg.add_task(calc_sum, kwargs=add_task_in_map.outputs.result)

# Set the final output of the workgraph
wg.outputs.result = sum_task.outputs.result

wg.run()

print(f"State of WorkGraph: {wg.state}")
print(f"Result: {wg.outputs.result.value}")

assert wg.outputs.result.value == 10

# Visualize the workgraph
wg.to_html()

# %%
# Conclusion
# ==========
# This tutorial has demonstrated how to construct workflows using the **node-graph programming paradigm** in `aiida-workgraph`.
# It presents an alternative approach to the Pythonic approach used in the rest of the documentation.
# The Pythonic approach serves as syntactic sugar that simplifies workflow construction, while node-graph programming  is the low-level approach where you build the graph piece by piece.
# This method offers maximum control but is more verbose and is generally reserved for advanced use cases, like programmatically generating a graph's structure.
