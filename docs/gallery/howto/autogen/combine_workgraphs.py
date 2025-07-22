"""
Combine workgraphs
==================
"""


# %%
# Introduction
# ------------
#
# In designing a complex workflow, it is often desired to reuse smaller, reusable components.
# In the following sections, we'll create a simple workgraph and integrate it into another.
#
# To start, let's define a workgraph to add two numbers and multiply the sum by a third.

from aiida_workgraph import WorkGraph, task
from aiida import load_profile

load_profile()

# %%
with WorkGraph("AddMultiply") as wg1:
    wg1.inputs = dict.fromkeys(["x", "y", "z"])
    wg1.outputs.result = (wg1.inputs.x + wg1.inputs.y) * wg1.inputs.z

wg1.to_html()

# %%
# We can see our two tasks, the linking of the sum to the first multiplication factor, and the assignment of the product as the final workgraph result.
# We're now ready to integrate our new **AddMultiply** workgraph into other workgraphs.

# %%
# Add a workgraph as a task
# -------------------------
#
# Adding a workgraph as a task of another is straightforward. We define a new workgraph, **AddMultiplyWorkGraph**, with a new task to generate a random number.
# We then call the **AddMultiply** workgraph as a task within this new workgraph, assigning the random number to the ``z`` input (the multiplication factor).
# Finally, we set the output of the **AddMultiply** workgraph as the output of the **AddMultiplyWorkGraph**.


@task
def generate_random_number(minimum, maximum):
    import random

    return random.randint(minimum, maximum)


with WorkGraph("AddMultiplyComposed") as wg2:
    wg2.inputs = dict.fromkeys(["min", "max", "x", "y"])

    wg2.outputs.result = wg1(
        inputs={
            "x": wg2.inputs.x,
            "y": wg2.inputs.y,
            "z": generate_random_number(
                minimum=wg2.inputs.min,
                maximum=wg2.inputs.max,
            ).result,
        }
    ).result

wg2.to_html()

# %%
# See how we're using **AddMultiply** as a regular task? It's as simple as that!
# Let's run our new workgraph and have a look at its result.

wg2.submit(
    inputs={
        "min": 1,
        "max": 10,
        "x": 1,
        "y": 2,
    },
    wait=True,
)

# %%
random_number_data_node = wg2.tasks.generate_random_number.outputs.result.value
final_result_data_node = wg2.outputs.result.value

print(f"Randomly generated number: {random_number_data_node.value}")
print(
    f"Final result: {final_result_data_node.value} = (1 + 2) * {random_number_data_node.value}"
)

# %%
# Let's have a look at the provenance graph.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg2.pk)

# %%
# Making it reusable
# -------------------
#
# Now, how do we make this reusable? Simple! We just wrap **AddMultiply** in a generating function.


def generate_add_multiply_workgraph():
    with WorkGraph("AddMultiply") as wg:
        wg.inputs = dict.fromkeys(["x", "y", "z"])
        wg.outputs.result = (wg.inputs.x + wg.inputs.y) * wg.inputs.z
    return wg


# %%
# Now we can generate an **AddMultiply** workgraph anytime we need one in our workflows.

# %%
# Summary
# -------
#
# Combining workgraphs is a straightforward process that allows for a modular approach to workflow design.
# By defining reusable components, we can easily integrate them into larger workflows, enhancing maintainability and readability.
# The example provided demonstrates how to create a simple workgraph and incorporate it into another, showcasing the flexibility of the AiiDA WorkGraph framework.
