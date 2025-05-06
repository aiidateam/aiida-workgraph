"""
WorkGraph
=================

This :class:`~aiida_workgraph.workgraph.WorkGraph` object is a collection of tasks and links.

Create and launch workgraph
----------------------------
"""

# %%
# First, create an empty workgraph:
#

from aiida_workgraph import WorkGraph, task
from aiida import load_profile

load_profile()

wg = WorkGraph(name="my_first_workgraph")

# %%
# Define and use tasks
#

# Define a task:
@task()
def add(x, y):
    return x + y


# Add tasks to the workgraph
add1 = wg.add_task(add, name="add1")
add2 = wg.add_task(add, name="add2")

# %%
# Add a link between tasks:

wg.add_link(add1.outputs.result, add2.inputs.x)
wg.to_html()

# %%
# Submit the workgraph:

wg.submit(inputs={"add1": {"x": 1, "y": 2}, "add2": {"y": 3}}, wait=True)

# %%
# Load workgraph from the AiiDA process
# -------------------------------------
#
# WorkGraph saves its data as an extra attribute in its process, allowing reconstruction of the WorkGraph from the process.

from aiida_workgraph import WorkGraph

wg_loaded = WorkGraph.load(wg.pk)

# %%
# Execute order
# -------------
# The tasks will be executed under the following conditions:
#
# - No input task
# - All input tasks are finished.
#
#
# Grouping Inputs and Outputs in a WorkGraph
# ------------------------------------------
# Defining **group-level** inputs and outputs allows you to:
#
# - Reuse inputs across multiple tasks (e.g., when several tasks share the same parameter).
# - Present only the necessary inputs to users, simplify the external interface of a complex workflow.
# - Collect and optionally rename outputs from individual tasks as grouped outputs.

wg = WorkGraph("test_workgraph_group_outputs")

# Define group-level input
wg.group_inputs.x = 2

# Add tasks using the group-level input
wg.add_task(add, "add1", x=wg.group_inputs.x, y=3)
wg.add_task(add, "add2", x=wg.group_inputs.x, y=wg.tasks.add1.outputs.result)

# Define group-level outputs to expose selected task results
wg.group_outputs.sum1 = wg.tasks.add1.outputs.result
wg.group_outputs.sum2 = wg.tasks.add2.outputs.result

# Run the workgraph
wg.submit(wait=True)

# Verify the final output
assert wg.group_outputs.sum2.value == 2 + (2 + 3)

# %%
# List of all Methods
# ----------------------------
#
# .. autoclass:: aiida_workgraph.workgraph.WorkGraph
#    :members:
