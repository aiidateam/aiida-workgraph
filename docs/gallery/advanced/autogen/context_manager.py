"""
==================================================
Write workflows using the context manager paradigm
==================================================
"""

# %%
# This guide introduces the **context manager paradigm** in `aiida-workgraph`.
# The paradigm provides the means to explicitly define a workflow as a linear, non-nested sequence of tasks.
# This includes conditional and iterative constructs, allowing the user to clearly see the logical flow of the workflow.
# By using this approach to writing workflows, the user gains the power to control tasks not yet executed, as all possible paths are fully displayed.
#
# In the following sections, we'll explore how to build workflows using the context manager paradigm.
# We'll cover simple sequential workflows, conditional branches, and iterative execution in loops.
# Along the way, we'll highlight the contrasts with the ``@task``-based approach as they come up.
#
# Let's get started!

# %%
# Creating a simple workflow
# ==========================
#
# We'll start by creating a simple arithmetic workflow using the context manager paradigm.
# We first recall how this is done using the ``@task`` decorator:

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


@task.graph
def AddMultiply(x, y, z):
    the_sum = add(x, y).result
    return multiply(the_sum, z).result


# %%
# Now, let's see how we can achieve the same result using the context manager paradigm.

from aiida_workgraph import WorkGraph


with WorkGraph("AddMultiplyContextManager") as wg:
    wg.inputs = dict.fromkeys(["x", "y", "z"])

    the_sum = add(
        x=wg.inputs.x,
        y=wg.inputs.y,
    ).result

    the_product = multiply(
        x=the_sum,
        y=wg.inputs.z,
    ).result

    wg.outputs.result = the_product


# %%
# A few things to note:
#
# - We explicitly name our workflow using the ``WorkGraph`` class
# - We get access to the ``WorkGraph`` instance (``wg``) directly
# - We use the ``WorkGraph`` context manager to define the workflow
# - Unlike the ``@task`` decorator, we explicitly set the inputs of the workflow using ``wg.inputs``
# - We pass to each task the inputs directly from the workflow's inputs
# - The results of the tasks are assigned to the workflow's outputs using ``wg.outputs.<socket_name>``
#
# How do they compare visually?

AddMultiply.build_graph(x=1, y=2, z=3).to_html()

# %%

wg.to_html()

# %%
# With the context manager approach, we now gain a clear representation of workflow inputs and their connections to tasks.
# This is one benefit of explicitly defining workflow inputs and assigning input sockets to tasks explicitly.
#
# .. admonition:: Take home message
#
#    Using the context manager paradigm allows for clear visualization of workflow inputs and their connections to tasks.


# %%
# Workflow inputs/outputs
# =======================
#
# As we saw in the previous section, we get more control over the workflow inputs when using the context manager paradigm.
# This is provided via direct access to the ``WorkGraph`` instance.
# Let's see how we can use this further.

# %%
# Nested namespaces
# -----------------
#
# We can define our workflow inputs (and outputs) using namespaces for clarity and convenience.
# For example, consider the following workflow:

with WorkGraph("AddThreeMultiplyContextManager") as wg:
    wg.inputs = {
        "add": {
            "first": dict.fromkeys(["x", "y"]),
            "second": dict.fromkeys(["x", "y"]),
        },
        "multiply": dict.fromkeys(["factor"]),
    }

    first_sum = add(
        x=wg.inputs.add.first.x,
        y=wg.inputs.add.first.y,
    ).result

    second_sum = add(
        x=wg.inputs.add.second.x,
        y=wg.inputs.add.second.y,
    ).result

    third_sum = add(
        x=first_sum,
        y=second_sum,
    ).result

    product = multiply(
        x=third_sum,
        y=wg.inputs.multiply.factor,
    ).result

    wg.outputs = {
        "sums": {
            "first": first_sum,
            "second": second_sum,
            "third": third_sum,
        },
        "product": product,
    }

# %%
# When defining inputs under namespaces (here, ``add`` and ``multiply``), we can access them using the dot notation when assigning to tasks.
# Similarly, we can access the outputs using the same notation.
# Let's run our workflow, now with a clear input layout:

from aiida import load_profile

load_profile()

wg.run(
    inputs={
        "add": {
            "first": {
                "x": 1,
                "y": 2,
            },
            "second": {
                "x": 3,
                "y": 4,
            },
        },
        "multiply": {
            "factor": 5,
        },
    },
)

print("\nResults:")
print("  Sums:")
print("    First:", wg.outputs.sums.first.value)
print("    Second:", wg.outputs.sums.second.value)
print("    Third:", wg.outputs.sums.third.value)
print("  Product:", wg.outputs.product.value)

# %%
# Let's see what this looks like visually.

wg.to_html()

# %%
# Again, we see our workflow inputs, this time clearly organized under the ``add`` and ``multiply`` namespaces, with clear connections to each task.
# The same applies to the outputs, which are grouped under the ``sums`` and ``product`` namespaces.

# %%
# .. admonition:: Take home message
#
#    Input/output namespaces provide a clear and convenient way to organize and access workflow data.

# %%
# Input metadata
# --------------
#
# Workflow inputs can also be defined using ``wg.add_input(...)``.
# The method allows you to provide an identifier (e.g. ``workgraph.int``) to the input, which is used for validation.

with WorkGraph() as wg:
    wg.add_input("workgraph.int", "x")  # validated as an integer
    wg.add_input("workgraph.int", "y")
    wg.add_input("workgraph.int", "z")
    ...

# %%
# .. tip::
#
#    When using the AiiDA GUI, providing an ``identifier`` to input sockets will associate the input with a GUI component, allowing users to interact with the input in a more user-friendly type-specific way (see this :ref:`GUI section <web-ui:detailed-socket-view>`).
#
# In the future, further details may be added when defining inputs, e.g., default values, descriptions, help messages, etc.

# %%
# Nested workflows
# ================
#
# We can reuse existing workflows as tasks within our workflow.
# Let's have a look at how this works in practice:


@task
def generate_random_number(minimum, maximum):
    import random

    return random.randint(minimum, maximum)


def generate_add_multiply_workgraph():
    with WorkGraph() as wg:
        wg.inputs = dict.fromkeys(["x", "y", "z"])

        the_sum = add(
            x=wg.inputs.x,
            y=wg.inputs.y,
        ).result

        the_product = multiply(
            x=the_sum,
            y=wg.inputs.z,
        ).result

        wg.outputs.result = the_product
    return wg


with WorkGraph("AddMultiplyComposed") as wg:
    wg.inputs = dict.fromkeys(["min", "max", "x", "y"])

    random_number = generate_random_number(
        minimum=wg.inputs.min,
        maximum=wg.inputs.max,
    ).result

    nested_wg = generate_add_multiply_workgraph()

    wg.outputs.result = nested_wg(
        inputs={
            "x": wg.inputs.x,
            "y": wg.inputs.y,
            "z": random_number,
        }
    ).result

# %%
# We define a new workgraph, `AddMultiplyComposed` that reuses an `AddMultiply` workgraph (here wrapped in a reusable generator function) with a random number as input.
# In our new workflow, we first call the generator function to get an instance of the workgraph.
# We then call it with its inputs (similar to ``.run(inputs=...)``).
# As a task, it returns a socket namespace, in which we defined a ``result`` socket.
# Finally, we assign this ``result`` socket as the ``result`` socket of our composed workflow.
#
# Let's have a look at the graph:

wg.to_html()

# %%
# You can compare this graph against the one generated by the decorator paradigm (see :ref:`here <howto:combine-workgraphs:add-workgraph>`).
# Note the presence of the workflow inputs and their connections to the various tasks (as discussed earlier).

# %%
# .. _advanced:context-manager:continue-workflow:
#
# Continuing a workflow
# =====================
#
# One of the key features of ``WorkGraph`` is its ability to continue previous workflows.
# When a workgraph finishes its execution, it saves its state in the AiiDA process node.
# This allows you to rebuild the workgraph from the process and add new tasks to continue the workflow.

# %%
# Continue with modified inputs
# -----------------------------
#
# Let's run an add-multiply workflow with a hardcoded multiplication factor:

with WorkGraph("AddMultiplyToBeContinued") as wg1:
    wg1.inputs = dict.fromkeys(["x", "y"])

    the_sum = add(
        x=wg1.inputs.x,
        y=wg1.inputs.y,
    ).result

    the_product = multiply(
        x=the_sum,
        y=3,
    ).result

    wg1.outputs = {
        "sum": the_sum,
        "product": the_product,
    }

wg1.run(
    inputs={
        "x": 1,
        "y": 2,
    },
)

print("\nResults:")
print(f"  Sum: {wg1.outputs.sum.value}")
print(f"  Product: {wg1.outputs.product.value}")

# %%
# Suppose we now want to rerun it with a different multiplication factor.
# Let's see how that's done:

with WorkGraph.load(wg1.pk) as wg2:
    wg2.name = "AddMultiplyModified"
    wg2.restart()
    wg2.tasks.multiply.inputs.y = 4

wg2.submit(wait=True)

print("\nResults:")
print(f"  Sum: {wg2.outputs.sum.value}")
print(f"  Product: {wg2.outputs.product.value}")

# %%
# Note that the sum has not changed (the ``value``, but more importantly, the ``pk``, as it is the same node).
# The product, however, is the result of the calculation repeating with the new input, hence a brand new node.

# %%
# Continue with new tasks
# -----------------------
#
# Let's now pick up the previous workgraph and extend it by a second addition, leveraging the results of the previous work.

with WorkGraph.load(wg2.pk) as wg3:
    wg3.name = "AddMultiplyContinued"
    wg3.inputs = dict.fromkeys(["z"])
    wg3.restart()
    new_sum = add(
        x=wg3.tasks.multiply.outputs.result,
        y=wg3.inputs.z,
    ).result
    wg3.outputs.new_sum = new_sum

wg3.to_html()

# %%
print(f"State of WorkGraph : {wg3.state}")
print(f"State of add       : {wg3.tasks.add.state}")
print(f"State of multiply  : {wg3.tasks.multiply.state}")
print(f"State of new add   : {wg3.tasks.add1.state}")

# %%
# Note the ``PLANNED`` new addition task. Let's run it.
# Let's run it with the new input:

wg3.submit(
    inputs={
        "z": 5,
    },
    wait=True,
)

print("\nResults:")
print(f"  Sum: {wg3.outputs.sum.value}")
print(f"  Product: {wg3.outputs.product.value}")
print(f"  New sum: {wg3.outputs.new_sum.value}")

# %%
# Again, note that the previous data nodes are the same.
# Only the new addition task ran and created a new data node.
#
# The full provenance remains intact, including the original workflow:

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg3.pk)

# %%
# Conclusion
# ==========
#
# In this section, you learned how to:
#
# - Use ``wg.inputs = {...}`` to define many inputs at once, or to define nested (namespaced) inputs to group related parameters
# - Use ``wg.add_input(...)`` to define a graph-level input and provide additional metadata (e.g., type validation)
# - Use graph-level inputs in tasks (``wg.inputs.<name>``)
# - Use ``wg.outputs.<name>`` to expose graph-level outputs from internal tasks
# - Load an existing workgraph using ``WorkGraph.load(pk)``
# - Continue a workgraph by modifying inputs or adding new tasks
