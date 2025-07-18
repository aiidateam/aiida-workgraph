"""
Continue a finished workgraph
=============================
"""


# %%
# Introduction
# ------------
#
# One of the key features of ``WorkGraph`` is its ability to continue previous workflows.
# When a workgraph finishes its execution, it saves its state in the AiiDA process node.
# This allows you to rebuild the workgraph from the process and add new tasks to continue the workflow.
#
# In the following sections, we will demonstrate how you can load a finished workgraph, amend it with new tasks, and continue the workflow.

from aiida_workgraph import WorkGraph
from aiida import load_profile

load_profile()

# %%
# Submit a workgraph
# ------------------
#
# We first create and submit a simply workgraph, so that we can reload and continue it later.

with WorkGraph("AddMultiplyToBeContinued") as wg1:
    wg1.inputs = dict.fromkeys(["x", "y"])
    the_sum = wg1.inputs.x + wg1.inputs.y
    the_product = the_sum * 3
    wg1.outputs.result = {
        "sum": the_sum,
        "product": the_product,
    }

wg1.submit(
    inputs={
        "x": 1,
        "y": 2,
    },
    wait=True,
)

# %%
print("Results:")
print(f"  Sum: {wg1.outputs.result.sum.value}")
print(f"  Product: {wg1.outputs.result.product.value}")

# %%
# Modify the workgraph and resubmit
# ---------------------------------
#
# It is possible to reload the workgraph, modify an input, and resubmit.
# We use the ``load`` method to load the previous workgraph by its ``pk``, mark it for restart, and modify the input of the multiplication task.
# This allows us to change the behavior of the workgraph without having to redefine it from scratch.

with WorkGraph.load(wg1.pk) as wg2:
    wg2.name = "AddMultiplyModified"
    wg2.restart()
    wg2.tasks.op_mul.inputs.y = 4

wg2.to_html()

# %%
wg2.submit(wait=True)

# %%
print("Results:")
print(f"  Sum: {wg2.outputs.result.sum.value}")
print(f"  Product: {wg2.outputs.result.product.value}")

# %%
# Note that the sum has not changed (the ``value``, but more importantly, the ``pk``, as it is the same node).
# The product, however, is the result of the calculation repeating with the new input, hence a brand new node.

# %%
# Continue a workgraph
# --------------------
#
# Let's now pick up the previous workgraph and extend it by a second addition, leveraging the results of the previous work.

with WorkGraph.load(wg2.pk) as wg3:
    wg3.name = "AddMultiplyContinued"
    wg3.inputs = dict.fromkeys(["z"])
    wg3.restart()
    wg3.outputs.result.new_sum = wg3.tasks.op_mul.outputs.result + wg3.inputs.z

wg3.to_html()

# %%
print(f"State of WorkGraph : {wg3.state}")
print(f"State of add       : {wg3.tasks.op_add.state}")
print(f"State of multiply  : {wg3.tasks.op_mul.state}")
print(f"State of new add   : {wg3.tasks.op_add1.state}")

# %%
# Note the ``PLANNED`` new addition task. Let's run it.

# %%
wg3.submit(
    inputs={
        "z": 5,
    },
    wait=True,
)

# %%
print("Results:")
print(f"  Sum: {wg3.outputs.result.sum.value}")
print(f"  Product: {wg3.outputs.result.product.value}")
print(f"  New sum: {wg3.outputs.result.new_sum.value}")

# %%
# Again, note that the previous data nodes are the same.
# Only the new addition task ran and created a new data node.

# %%
# The extended provenance
# -----------------------
#
# The full provenance remains intact, including the original workflow.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg3.pk)

# %%
# Summary
# -------
#
# In this section, you learned how to load an existing workflow and rerun it with modified input, or continue it with new tasks.
# These features allow you to extend workflows without losing the original provenance, enabling a flexible and iterative approach to workflow management in AiiDA.
