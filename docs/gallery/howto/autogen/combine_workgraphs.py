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

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


@task.graph
def AddMultiply(x, y, z):
    the_sum = add(x=x, y=y).result
    return multiply(x=the_sum, y=z).result


# %%
# We're now ready to integrate our new `AddMultiply` workgraph into other workgraphs.

# %%
# .. _howto:combine-workgraphs:add-workgraph:
#
# Add a workgraph as a task
# -------------------------
#
# Adding a workgraph as a task of another is straightforward. We define a new workgraph, `AddMultiplyComposed`, with a new task to generate a random number.
# We then call the `AddMultiply` workgraph as a task within this new workgraph, assigning the random number to the ``z`` input (the multiplication factor).
# We return the output of the `AddMultiply` workgraph, effectively assigning it as the output of the `AddMultiplyComposed`.


@task
def generate_random_number(minimum, maximum):
    import random

    return random.randint(minimum, maximum)


@task.graph
def AddMultiplyComposed(minimum, maximum, x, y):
    random_number = generate_random_number(minimum=minimum, maximum=maximum).result
    return AddMultiply(x=x, y=y, z=random_number).result


wg = AddMultiplyComposed.build(minimum=1, maximum=10, x=1, y=2)

# %%
# See how we're using `AddMultiply` as a regular task? It's as simple as that!
# This is also clear in when we visualize the workgraph:

wg.to_html()

# %%
# .. tip::
#
#    When using the AiiDA GUI, you can inspect a nested workgraph by clicking on the task node, then on **Go to WorkGraph**.
#    To learn more about this, see this :ref:`GUI section<web-ui:nested-workgraphs>`.
#
# Let's run our composed workgraph and have a look at its result:

from aiida import load_profile

load_profile()

wg.run()

random_number = wg.tasks.generate_random_number.outputs.result.value.value

print("\nResults:")
print("  Random number:", random_number)
print("  Final result:", f"{wg.outputs.result.value.value} = (1 + 2) * {random_number}")

# %%
# Let's have a look at the provenance graph:


wg.generate_provenance_graph()

# %%
# Summary
# -------
#
# Combining workgraphs is a straightforward process that allows for a modular approach to workflow design.
# By defining reusable components, we can easily integrate them into larger workflows, enhancing maintainability and readability.
# The example provided demonstrates how to create a simple workgraph and incorporate it into another, showcasing the flexibility of the AiiDA WorkGraph framework.
