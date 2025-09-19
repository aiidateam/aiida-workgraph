"""
=================
WorkGraph
=================

WorkGraph is a collection of tasks and links.


"""

# %%
# Create workgraph
# ========================
# First, create an empty workgraph:
#

from aiida_workgraph import WorkGraph, task

wg = WorkGraph(name='my_first_workgraph')

# %%
# Define and use tasks
#


# Define a task:
@task()
def add(x, y):
    return x + y


# Add tasks to the workgraph
add1 = wg.add_task(add, name='add1')
add2 = wg.add_task(add, name='add2')

# %%
# Add a link between tasks:

wg.add_link(add1.outputs.result, add2.inputs.x)

# Visualize the graph
wg.to_html()

# %%
# Execute the workgraph
# ========================
# With the graph defined, you can now execute it. You provide the inputs for the tasks.

from aiida import load_profile

load_profile()
wg.run(inputs={'add1': {'x': 1, 'y': 2}, 'add2': {'y': 3}})

# %%
# Graph-level inputs and outputs
# ========================================
# As workflows grow, managing inputs for many tasks can become cumbersome.
# WorkGraph allows you to define **graph-level** inputs and outputs to create a
# cleaner, more user-friendly interface for your complex logic.
#
# This lets you:
#
# - **Reuse** a single input across multiple tasks.
# - **Hide** internal complexity and only expose essential inputs.
# - **Collect** and rename key results as named workflow outputs.

wg = WorkGraph('graph_inputs_outputs')

# Define graph-level input
wg.inputs.x = 2

# Add tasks using the graph-level input
wg.add_task(add, 'add1', x=wg.inputs.x, y=3)
wg.add_task(add, 'add2', x=wg.inputs.x, y=wg.tasks.add1.outputs.result)

# Define graph-level outputs to expose selected task results
wg.outputs.sum1 = wg.tasks.add1.outputs.result
wg.outputs.sum2 = wg.tasks.add2.outputs.result

# Run the WorkGraph
wg.run()

# Verify the final output
assert wg.outputs.sum2.value == 2 + (2 + 3)

# Visualize the graph with inputs and outputs
wg.to_html()

# %%
# Context variables
# ========================================
# Context variables (`ctx`) are used to store and pass intermediate data within a
# workflow that isn't directly an input or output of a task. This is
# especially useful for workflows with conditional logic (if/else) or loops,
# where you need to manage state between steps.

wg = WorkGraph(name='context_example')
# Setting the ``ctx`` attribute of the WorkGraph directly, on initialization
wg.ctx = {'x': 2, 'data.y': 3}
wg.add_task(add, 'add1', x=wg.ctx.x, y=wg.ctx.data.y)
# Assign the result of a task to a context variable
wg.ctx.sum = wg.tasks.add1.outputs.result
# Use the context variable in another task
wg.add_task(add, 'add2', x=wg.ctx.x, y=wg.ctx.sum)

# %%
# Context variables can be nested, allowing you to organize complex data structures.
# For example, you can store multiple results in a structured way:

wg.ctx.data = {
    'sum1': wg.tasks.add1.outputs.result,
    'sum2': wg.tasks.add2.outputs.result,
}

# %%
# WorkGraph engine
# ========================================
# The WorkGraph engine operates on a **dataflow programming** model. Once
# submitted, the engine continuously monitors the tasks in the graph. A task
# is executed only when all of its inputs are available. This means:
#
# - Tasks with no inputs are executed first.
# - A task starts only after all the upstream tasks linked to its inputs are finished.
#
