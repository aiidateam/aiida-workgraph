"""
========================================
Sockets
========================================

In WorkGraph, a **Socket** is a fundamental concept that defines the connection points for data to flow between tasks. Think of them as the inputs and outputs of each step in your workflow. They are analogous to AiiDA's `Ports` but with an important extension: sockets are designed to manage not just the data type but also the **links** between tasks in the graph.

This guide will walk you through how to define, customize, and organize sockets for your tasks.

"""

from aiida.manage import load_profile
from aiida_workgraph import task, WorkGraph, spec
from aiida import orm
from typing import Any

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
print('Input sockets: ', task1.get_input_names())
print('Output sockets: ', task1.get_output_names())


# %%
# Customizing Output Sockets
# ==============================
# When a function needs to return multiple outputs, you can define custom output sockets with the outputs argument in the ``@task.decorator``. The task will map returned values to the named sockets in two ways:
#
# - The function returns a dictionary: The keys of the returned dictionary must match the socket names.


@task
def add_and_subtract(x, y) -> spec.namespace(sum=Any, difference=Any):
    """Return the sum and difference of two numbers in a dict."""
    return {
        'sum': x + y,
        'difference': x - y,
    }


# Inspect the new input and output sockets
task2 = wg.add_task(add_and_subtract, x=3, y=4)
print('Input sockets: ', task2.get_input_names())
print('Output sockets: ', task2.get_output_names())

# %%
# - The function returns a tuple: The elements of the tuple are mapped to the sockets in the order they are declared in the outputs definition.
#
# .. note::
#
#    Be sure that the number of elements in the returned tuple matches the number of defined output sockets.


@task
def add_and_subtract(x, y) -> spec.namespace(sum=Any, difference=Any):
    """Return the sum and difference of two numbers as a tuple."""
    return x + y, x - y


# Inspect the new input and output sockets
task2 = wg.add_task(add_and_subtract, x=3, y=4)
print('Input sockets: ', task2.get_input_names())
print('Output sockets: ', task2.get_output_names())


# %%
# The ``identifier`` is used to indicates the data type. The data
# type tell the code how to display the port in the GUI, validate the data,
# and serialize data into database.
# We use ``workgraph.Any`` for any data type.
#
# .. warning::
#   **The data validation feature is experimental.** The API is subject to change in future releases. We welcome your feedback on its functionality.
#
# For convenience, we also support using a list of strings as the output definition, which is equivalent to using ``workgraph.Any`` for each output socket.


@task(outputs=['sum', 'difference'])
def add_and_subtract(x, y):
    """This task returns both the sum and difference of two numbers."""
    return {'sum': x + y, 'difference': x - y}


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
print('Input x: ', task3.inputs.x)


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
    print(f'3 to the default power of 2 is: {outputs1.result.value}')
    print(f'3 to the power of 3 is: {outputs2.result.value}')


# %%
# Organizing Sockets with Namespaces
# =====================================
# As workflows grow, you might have many related inputs or outputs. To keep them organized and avoid name clashes, you can group them into a **namespace**.
#
# Simple Output Namespace
# -----------------------
# If a task has multiple outputs (as defined in the `outputs` decorator argument), WorkGraph automatically groups them into a `SocketNamespace`.

with WorkGraph('simple_namespace_example') as wg:
    # The `add_and_subtract` task from before has two outputs.
    # The `outputs` object becomes a namespace.
    outputs = add_and_subtract(x=10, y=4)
    wg.run()

    # You can access the results like attributes of an object:
    print('Accessing outputs from a namespace:')
    print(f'  Sum: {outputs.sum.value}')
    print(f'  Difference: {outputs.difference.value}')


# %%
# Nested Namespaces
# -----------------
# For more complex data structures, you can define nested namespaces. This allows you to create a hierarchical organization for your sockets.

out = spec.namespace(
    normal=spec.namespace(sum=Any, product=Any),
    squared=spec.namespace(sum=Any, product=Any),
)


@task(outputs=out)
def advanced_math(x, y):
    """A task with a nested output structure."""
    return {
        'normal': {'sum': x + y, 'product': x * y},
        'squared': {'sum': x**2 + y**2, 'product': x**2 * y**2},
    }


with WorkGraph('nested_namespace_example') as wg:
    outputs = advanced_math(x=2, y=3)
    wg.run()
    print('\nAccessing outputs from a nested namespace:')
    print(f'  Normal sum: {outputs.normal.sum.value}')
    print(f'  Squared product: {outputs.squared.product.value}')
    assert outputs.normal.sum.value == 5
    assert outputs.squared.product.value == 36


# %%
# .. _dynamic_namespaces:
# Dynamic namespaces
# -----------------------------------------------------------
# In many workflows, you might not know the exact number of outputs a task
# will generate until it runs. A **dynamic namespace** is a powerful feature
# in AiiDA-WorkGraph that allows a task to emit an arbitrary number of named
# output ports at runtime.
#
# To create one, define a namespace socket and set its metadata to `{"dynamic": True}`.


@task
def generate_squares(n: int) -> spec.namespace(squares=spec.dynamic(int)):
    """Generates a dynamic number of square values."""
    # The return dictionary must match the output socket names.
    # The value for the "squares" dynamic namespace is another dictionary,
    # where each key-value pair will become a separate AiiDA node.
    return {'squares': {f'n_{i}': i**2 for i in range(n)}}


# %%
# This approach offers two significant advantages:
#
# 1.  **Fine-grained provenance**: Each key-value pair (`"n_0": 0`, `"n_1": 1`, etc.)
#     is stored as an individual AiiDA node. This gives you a complete and transparent
#     graph, allowing you to trace the history and reuse of every single result.
#
# 2.  **Database integrity**: It avoids serializing large, complex Python objects
#     into a single, opaque file.
#     Instead, each result is a native, queryable entry in the AiiDA database.
#
# To consume these dynamic outputs, we can create a task that accepts an arbitrary
# number of keyword arguments.


@task
def sum_all(**data):
    """A task to sum all values from a dictionary of inputs."""
    # The **data syntax gathers all keyword arguments into a dictionary named 'data'.
    return sum(data.values())


# %%
# Now, let's build the WorkGraph.
with WorkGraph('dynamic_namespace_example') as wg:
    # The first task generates a dynamic set of outputs under the "squares" namespace.
    dynamic_outputs = generate_squares(n=4)

    # The entire "squares" namespace is then linked to the next task.
    total = sum_all(data=dynamic_outputs.squares)

    wg.run()

    # You can access the individual dynamic outputs
    print('\nIndividual dynamic outputs:')
    for i in range(4):
        print(f'  n_{i}: {dynamic_outputs.squares[f"n_{i}"].value}')

    print(f'\nSum of all dynamic outputs: {total.result.value}')
    assert total.result.value == 14

# %%
# .. important::
#
#    **Connecting sockets: why `data=...` is essential**
#
#    When linking a dynamic namespace to a task that accepts keyword arguments (`**kwargs`),
#    you **must** explicitly name the input argument. In our example:
#
#    .. code-block:: python
#
#        # CORRECT: Connect the 'squares' output socket to the 'data' input socket.
#        sum_all(data=dynamic_outputs.squares)
#
#    This is because you are not passing values directly, but rather connecting the *sockets*
#    of the graph. `dynamic_outputs.squares` is a reference to the *output socket* that will
#    hold all the dynamic results. The task `sum_all` has a corresponding *input socket*
#    named `data` (derived from the `**data` signature). The line above tells WorkGraph to
#    wire the entire dynamic namespace from the output socket to this specific input socket.
#
#    Using Python's argument unpacking syntax will **fail**:
#
#    .. code-block:: python
#
#        # INCORRECT: This will raise an error.
#        sum_all(**dynamic_outputs.squares)
#
#    At the time the graph is built, `dynamic_outputs.squares` is a socket object, not a
#    Python dictionary with concrete values. The `**` operator attempts to unpack a dictionary
#    immediately, which cannot be done on a socket reference that represents future results.
#    You must name the argument to establish the connection for the engine to resolve later.
