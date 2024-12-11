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

wg = WorkGraph(name="my_first_workgraph")

# %%
# Define and use tasks
#

# Define a task using a `calcfunction`:
@task.calcfunction()
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
# Group outputs
# -------------
# You can output the results of the tasks as the output of the WorkGraph.

wg = WorkGraph("test_workgraph_group_outputs")
wg.add_task(add, "add1", x=2, y=3)
wg.group_outputs = [{"name": "sum", "from": "add1.result"}]
wg.submit(wait=True)
assert wg.process.outputs.sum.value == 5

# %%
# List of all Methods
# ----------------------------
#
# .. autoclass:: aiida_workgraph.workgraph.WorkGraph
#    :members:
