"""
========================================
Sockets
========================================

In WorkGraph, a **Socket** is a fundamental concept that defines the connection points for data to flow between tasks. Think of them as the inputs and outputs of each step in your workflow. They are analogous to AiiDA's `Ports` but with an important extension: sockets are designed to manage not just the data type but also the **links** between tasks in the graph.

This guide will walk you through how to define, customize, and organize sockets for your tasks.

"""
from aiida.manage import load_profile
from aiida_workgraph import task, WorkGraph
from aiida import orm

# Load the AiiDA profile to interact with the database
load_profile()

# %%
# Automatic Socket Creation
# ===========================
# The easiest way to create sockets is to let WorkGraph do it for you. When you decorate a Python function as a ``@task``, WorkGraph automatically inspects its signature:
#
# * Function arguments become **input sockets**.
# * The return value becomes a single **output socket** named ``result``.
#
# Beside, there are also some built-in sockets from the WorkGraph, like ``_wait`` and ``_outputs``.


@task()
def multiply(x, y):
    """A simple task to multiply two numbers."""
    return x * y


# Let's inspect the automatically generated sockets
wg = WorkGraph()
task1 = wg.add_task(multiply, x=3, y=4)
print("Input sockets: ", task1.get_input_names())
print("Output sockets: ", task1.get_output_names())


# %%
# Customizing Output Sockets
# ==============================
# When a function needs to return multiple outputs, you can define custom output sockets with the outputs argument in the ``@task.decorator``. The task will map returned values to the named sockets in two ways:
#
# - The function returns a dictionary: The keys of the returned dictionary must match the socket names.


@task(
    outputs={
        "sum": {"identifier": "workgraph.Any"},
        "difference": {"identifier": "workgraph.Any"},
    }
)
def add_and_subtract(x, y):
    """Return the sum and difference of two numbers in a dict."""
    return {
        "sum": x + y,
        "difference": x - y,
    }


# Inspect the new input and output sockets
task2 = wg.add_task(add_and_subtract, x=3, y=4)
print("Input sockets: ", task2.get_input_names())
print("Output sockets: ", task2.get_output_names())

# %%
# - The function returns a tuple: The elements of the tuple are mapped to the sockets in the order they are declared in the outputs definition.
#
# .. note::
#
#    Be sure that the number of elements in the returned tuple matches the number of defined output sockets.


@task(
    outputs={
        "sum": {"identifier": "workgraph.Any"},
        "difference": {"identifier": "workgraph.Any"},
    }
)
def add_and_subtract(x, y):
    """Return the sum and difference of two numbers as a tuple."""
    return x + y, x - y


# Inspect the new input and output sockets
task2 = wg.add_task(add_and_subtract, x=3, y=4)
print("Input sockets: ", task2.get_input_names())
print("Output sockets: ", task2.get_output_names())


# %%
# The ``identifier`` is used to indicates the data type. The data
# type tell the code how to display the port in the GUI, validate the data,
# and serialize data into database.
# We use ``workgraph.Any`` for any data type. For the moment, the data validation is
# experimentally supported, and the GUI display is not implemented. Thus,
# I suggest you to always ``workgraph.Any`` for the port.
#
# For convenience, we also support using a list of strings as the output definition, which is equivalent to using ``workgraph.Any`` for each output socket.


@task(outputs=["sum", "difference"])
def add_and_subtract(x, y):
    """This task returns both the sum and difference of two numbers."""
    return {"sum": x + y, "difference": x - y}


# %%
# Socket Types from Python Type Hints
# =========================================
# To ensure data integrity and leverage AiiDA's data provenance, sockets have types. WorkGraph can automatically assign AiiDA-compatible socket types based on standard Python type hints in your function signature. This is the recommended modern approach.
#
# =================== ==================
#  Python Type Hint    AiiDA Socket Type
# =================== ==================
#  ``int``             ``orm.Int``
#  ``float``           ``orm.Float``
#  ``str``             ``orm.Str``
#  ``bool``            ``orm.Bool``
#  ``list``            ``orm.List``
#  ``dict``            ``orm.Dict``
# =================== ==================
#


@task.calcfunction()
def add_typed(x: orm.Int, y: orm.Float) -> orm.Float:
    """A task with typed sockets to add an Int and a Float."""
    return orm.Float(x.value + y.value)


# When you call this task, WorkGraph expects AiiDA data types
task3 = wg.add_task(add_typed, x=orm.Int(3), y=orm.Int(4))
print("Input x: ", task3.inputs.x)


# %%
# Default Values as Socket Properties
# -----------------------------------------
# A socket can have a **property**, which is a default value used when no data is linked to an input. The most Pythonic way to define a property is by using a default argument in your function definition.


@task
def power(base: int, exponent: int = 2):
    """Calculates base to the power of exponent, defaulting to square."""
    return base**exponent


# The `exponent` socket now has a default value of 2.
# We only need to provide the `base`.
with WorkGraph() as wg:
    outputs1 = power(base=3)
    # We can still override the default value by providing an input.
    outputs2 = power(base=3, exponent=3)
    wg.run()
    print(f"3 to the default power of 2 is: {outputs1.result.value}")
    print(f"3 to the power of 3 is: {outputs2.result.value}")


# %%
# Organizing Sockets with Namespaces
# =====================================
# As workflows grow, you might have many related inputs or outputs. To keep them organized and avoid name clashes, you can group them into a **namespace**.
#
# Simple Output Namespace
# -----------------------
# If a task has multiple outputs (as defined in the `outputs` decorator argument), WorkGraph automatically groups them into a `SocketNamespace`.

with WorkGraph("simple_namespace_example") as wg:
    # The `add_and_subtract` task from before has two outputs.
    # The `outputs` object becomes a namespace.
    outputs = add_and_subtract(x=10, y=4)
    wg.run()

    # You can access the results like attributes of an object:
    print(f"Accessing outputs from a namespace:")
    print(f"  Sum: {outputs.sum.value}")
    print(f"  Difference: {outputs.difference.value}")


# %%
# Nested Namespaces
# -----------------
# For more complex data structures, you can define nested namespaces. This allows you to create a hierarchical organization for your sockets.


@task(
    outputs={
        "normal": {
            "identifier": "workgraph.namespace",
            "sockets": {
                "sum": {"identifier": "workgraph.any"},
                "product": {"identifier": "workgraph.any"},
            },
        },
        "squared": {
            "identifier": "workgraph.namespace",
            "sockets": {
                "sum": {"identifier": "workgraph.any"},
                "product": {"identifier": "workgraph.any"},
            },
        },
    }
)
def advanced_math(x, y):
    """A task with a nested output structure."""
    return {
        "normal": {"sum": x + y, "product": x * y},
        "squared": {"sum": x**2 + y**2, "product": x**2 * y**2},
    }


with WorkGraph("nested_namespace_example") as wg:
    outputs = advanced_math(x=2, y=3)
    wg.run()
    print("\nAccessing outputs from a nested namespace:")
    print(f"  Normal sum: {outputs.normal.sum.value}")
    print(f"  Squared product: {outputs.squared.product.value}")
    assert outputs.normal.sum.value == 5
    assert outputs.squared.product.value == 36


# %%
# Dynamic Namespaces
# ------------------
# Sometimes, you don't know the number of outputs a task will generate beforehand. A **dynamic namespace** can accept a variable number of sockets at runtime.
#
# To create one, set `"metadata": {"dynamic": True}` on a namespace socket.


@task(
    outputs={
        "squares": {
            "identifier": "workgraph.namespace",
            "metadata": {"dynamic": True},
        }
    }
)
def generate_squares(n: int):
    """Generates a dynamic number of square values."""
    return {"squares": {f"n_{i}": i**2 for i in range(n)}}


@task
def sum_all(**inputs):
    """A task to sum all values from a dictionary."""
    return sum(inputs.values())


with WorkGraph("dynamic_namespace_example") as wg:
    # The first task generates a dynamic set of outputs under the "squares" namespace.
    dynamic_outputs = generate_squares(n=4)

    # The entire dynamic namespace `dynamic_outputs.squares` is linked as a
    # single dictionary input to the next task.
    total = sum_all(inputs=dynamic_outputs.squares)
    wg.run()

    # You can access the individual dynamic outputs
    print(f"\nIndividual dynamic outputs:")
    for i in range(4):
        print(f"  n_{i}: {dynamic_outputs.squares[f'n_{i}'].value}")

    print(f"Sum of all dynamic outputs: {total.result.value}")
    assert total.result.value == 14
