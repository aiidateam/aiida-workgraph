"""
==================================
Nested and dynamic namespace
==================================
"""
# %%
# Introduction
# ============
# In aiida-core, the input or output is represents by `port`.
# The namespaces can be arbitrarily nested with ports and so are called port namespaces.
# The port`namespace` is widely used because:
# - keep the interface readable as the number of inputs grows
# - compose workflows without name clashes
# - dynamic namespace can accept unknown number of inputs

# same for the outputs.
# In aiida-workgraph, we use Socket to wrap the port to let it also allow
# - support link between the inputs and outputs
# And similary, the `SocketNamespace` is used to represent the namespace port.


# %%
# Simple example of namespace output
# ==================================

# Load the AiiDA profile.
from aiida import load_profile
from aiida_workgraph import WorkGraph, task

load_profile()


@task(outputs=["sum", "product"])
def add_and_multiply(x, y):
    """Add and multiply two numbers."""
    return {"sum": x + y, "product": x * y}


with WorkGraph() as wg:
    # The `add_and_multiply` task is called with two inputs.
    # The outputs are automatically wrapped in a `SocketNamespace`.
    outputs = add_and_multiply(x=1, y=2)
    wg.run()
    print("outputs: ", outputs)


# %%
# Here the `outputs` of the task is a "SocketNamespace" which includes two sockets:
# - `sum` and `product`.
# One can access the outputs by the name of the socket:

print("sum:", outputs.sum)
print("product:", outputs.product)


# %%
# Nested namespace
# ==================
# The `SocketNamespace` can be nested to represent a more complex structure.


@task(
    outputs=[
        {
            "name": "normal",
            "identifier": "workgraph.namespace",
            "sockets": {
                "sum": {"identifier": "workgraph.any"},
                "product": {"identifier": "workgraph.any"},
            },
        },
        {
            "name": "square",
            "identifier": "workgraph.namespace",
            "sockets": {
                "sum": {"identifier": "workgraph.any"},
                "product": {"identifier": "workgraph.any"},
            },
        },
    ]
)
def add_and_multiply_nested(x, y):
    """Add and multiply two numbers."""
    return {
        "normal": {"sum": x + y, "product": x * y},
        "square": {"sum": x**2 + y**2, "product": x**2 * y**2},
    }


with WorkGraph() as wg:
    # The `add_and_multiply_nested` task is called with two inputs.
    # The outputs are automatically wrapped in a `SocketNamespace`.
    outputs = add_and_multiply_nested(x=1, y=2)
    wg.run()
    print("outputs: ", outputs)
    print("normal sum:", outputs.normal.sum)
    print("normal product:", outputs.normal.product)
    print("square sum:", outputs.square.sum)
    print("square product:", outputs.square.product)


# %%
# Dynamic namespace
# ==================
# The `SocketNamespace` can also be dynamic, which means that the number of inputs or outputs is not fixed.


@task()
def generate_square(n):
    """Generate a list of integers from 0 to n-1."""
    result = {f"square_{i}": i**2 for i in range(n)}
    return result


with WorkGraph() as wg:
    # The `generate_square` task is called with one input.
    # The outputs are automatically wrapped in a `SocketNamespace`.
    outputs = generate_square(n=5)
    wg.run()
    print("outputs: ", outputs)
    for i in range(5):
        print(f"square_{i}:", outputs.result[f"square_{i}"])  # access dynamic output
