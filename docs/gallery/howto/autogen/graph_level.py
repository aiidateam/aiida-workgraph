"""
Define graph-level inputs and outputs
=====================================
"""


# %%
# Introduction
# ------------
#
# When constructing complex workflows, you may encounter tasks that share input parameters (e.g. ``code``).
# Also, you may wish to reorder or rename internal task outputs at the top level (e.g. ``wg.outputs.optimized_stuff`` instead of ``wg.outputs.optimize.stuff``).
#
# ``WorkGraph`` allows you to define **graph-level inputs and outputs** to:
#
# - **Share inputs** across multiple tasks in a graph
# - **Aggregate outputs** from internal tasks, optionally organizing or renaming them
# - **Simplify the interface** by exposing shared parameters to users (retaining the flexibility of internal task parameter definitions)

from aiida_workgraph import WorkGraph, task
from aiida import load_profile

load_profile()

# %%
# Defining graph-level inputs and outputs
# ---------------------------------------

with WorkGraph("GraphLevelInput") as wg:
    # Define graph-level inputs
    wg.inputs = dict.fromkeys(["x", "y", "z"])

    # Define the tasks
    the_sum = wg.inputs.x + wg.inputs.y
    the_product = the_sum * wg.inputs.z

    # Define graph-level outputs
    wg.outputs.sum = the_sum
    wg.outputs.product = the_product

wg.to_html()

# %%
wg.submit(
    inputs={
        "x": 1,
        "y": 2,
        "z": 3,
    },
    wait=True,
)

# %%
print("\nGraph-level outputs:")
print("  Sum:", wg.outputs.sum.value)
print("  Product:", wg.outputs.product.value)

# %%
# Providing graph-level inputs metadata
# -------------------------------------
#
# Graph-level inputs can also be defined using ``wg.add_input(...)``.
# The method allows you to provide an identifier (e.g. ``workgraph.int``) to the input, which is used for validation.
# Also, if you use the AiiDA GUI, providing an ``identifier`` will associate the input with a GUI component, allowing users to interact with the input in a more user-friendly type-specific way.

with WorkGraph("GraphLevelInputMetadata") as wg:
    wg.add_input("workgraph.int", "x")  # validated as an integer
    wg.add_input("workgraph.int", "y")
    wg.add_input("workgraph.int", "z")
    wg.outputs.result = (wg.inputs.x + wg.inputs.y) * wg.inputs.z

# %%
# In the future, further details may be added when defining inputs, e.g., default values, descriptions, help messages, etc.

# %%
# Nested graph-level inputs and outputs
# -------------------------------------
#
# Graph-level inputs and outputs can be **nested** allowing you to group related parameters and results.
# Here we're using the ``add_task`` interface to more clearly define the names of our tasks.


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


with WorkGraph("GraphLevelInputOutputNested") as wg:
    wg.inputs = {
        "add": {
            "first": dict.fromkeys(["x", "y"]),
            "second": dict.fromkeys(["x", "y"]),
        },
        "multiply": dict.fromkeys(["factor"]),
    }

    first_sum = wg.inputs.add.first.x + wg.inputs.add.first.y
    second_sum = wg.inputs.add.second.x + wg.inputs.add.second.y
    third_sum = first_sum + second_sum
    product = third_sum * wg.inputs.multiply.factor

    wg.outputs.results = {
        "sums": {
            "first": first_sum,
            "second": second_sum,
            "third": third_sum,
        },
        "product": product,
    }


wg.to_html()

# %%
# We can now run our workgraph with a clear input layout.
#
# .. note:
#
#    ``WorkGraph`` will automatically serialize the raw Python data into the corresponding AiiDA Data nodes (e.g., an ``int`` becomes ``orm.Int``, a ``str`` becomes ``orm.Str``, etc.) before execution.
# The exact serialization logic and all supported types (and how to register your own custom serializers) are described in detail in the **Data Serialization** section.

wg.submit(
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
    wait=True,
)

# %%
print("\nResults:")
print("  Sums:")
print("    First:", wg.outputs.results.sums.first.value)
print("    Second:", wg.outputs.results.sums.second.value)
print("    Third:", wg.outputs.results.sums.third.value)
print("  Product:", wg.outputs.results.product.value)

# %%
# When we inspect the outputs of the workgraph process, we see ``sums`` and ``product`` are grouped under the ``results`` output.

import subprocess

subprocess.run(["verdi", "process", "show", str(wg.pk)], check=True)

# %%
# Summary
# -------
#
# In this section, you learned how to:
#
# - Use ``wg.inputs = {...}`` to define many inputs at once, or to define nested (namespaced) inputs to group related parameters
# - Use ``wg.add_input(...)`` to define a graph-level input and provide additional metadata (e.g., type validation)
# - Use graph-level inputs in tasks (``wg.inputs.<name>``)
# - Use ``wg.outputs.<name>`` to expose graph-level outputs from internal tasks
