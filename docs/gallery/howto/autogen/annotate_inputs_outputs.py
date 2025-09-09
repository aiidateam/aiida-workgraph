"""
Use annotations to control data provenance
==========================================
"""

# %%
# Introduction
# ------------
#
# In ``aiida-workgraph``, a critical feature is tracking task inputs and outputs to ensure data provenance and reproducibility.
# To achieve this, task functions must be annotated with input and output specifications.
# This tells the WorkGraph how to handle, serialize, and store data as individual AiiDA nodes.
#
# This process addresses two key aspects of provenance:
#
# - **Data creation**: How should data be created and stored? For example, if a task returns a nested dictionary, should it be stored as a single entity or unpacked into separate nodes?
# - **Data lineage**: Where does the data come from? How does it flow between tasks in the workflow?
#
# This guide will walk you through the various ways to annotate your tasks.

import typing as t
from aiida import load_profile
from aiida_workgraph import task
from aiida_workgraph.utils import generate_node_graph, get_process_summary

load_profile()


# %%
# Data creation
# -------------
#
# Data creation is controlled at the ``Calculation`` level.
# In AiiDA, a ``Calculation`` is a process that performs a computation and creates new data nodes.
#
# Static namespaces for outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Let's define two tasks that perform the same calculation but have different output annotations.
#
# The first task, ``add_multiply1``, has no output specification.
# Consequently, the returned dictionary will be stored as a single AiiDA ``Dict`` node.
# The second task, ``add_multiply2``, uses an output namespace to specify that each key-value pair in the dictionary should be stored as a separate AiiDA node.

# Import the `namespace` module for specifications
from aiida_workgraph import namespace


@task
def add_multiply1(x, y):
    """Return a dictionary, which will be stored as a single Dict node."""
    return {"sum": x + y, "product": x * y}


@task
def add_multiply2(
    x: int,
    y: int,
) -> t.Annotated[
    dict,
    namespace(sum=int, product=int),
]:
    """Return a dictionary, but its elements are stored as separate Int nodes."""
    return {"sum": x + y, "product": x * y}


@task.graph
def add_multiply_graph(x: int, y: int):
    """A graph to run both versions of the add_multiply task."""
    add_multiply1(x=x, y=y)
    add_multiply2(x=x, y=y)


wg = add_multiply_graph.build_graph(x=1, y=2)
wg.run()

# %%
#
# .. note::
#
#    In ``add_multiply2``, we also annotated the input types (``x: int, y: int``).
#    This adds a layer of validaton to ensure only integers are passed to the task.
#    This feature is experimental - its API and behavior may change in future releases.
#
# Let's visualize the data provenance of our executed workflow:

generate_node_graph(wg.pk)

# %%
# As the provenance graph shows, ``add_multiply1`` has a single output node (``result``), while ``add_multiply2`` has two separate output nodes (``sum`` and ``product``), as defined in its namespace.
#
# Static namespaces for inputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Similarly, we can annotate inputs to unpack a dictionary into distinct data nodes.
# We use Python's standard ``typing.Annotated`` to attach the `aiida-workgraph` namespace metadata to the input type.


@task
def add_multiply3(
    data: t.Annotated[
        dict,
        namespace(x=int, y=int),
    ],
) -> t.Annotated[
    dict,
    namespace(sum=int, product=int),
]:
    """Take a dictionary as input but treat 'x' and 'y' as separate nodes."""
    return {"sum": data["x"] + data["y"], "product": data["x"] * data["y"]}


@task.graph
def add_multiply_graph_inputs(x: int, y: int):
    return add_multiply3(data={"x": x, "y": y})


wg = add_multiply_graph_inputs.build_graph(x=1, y=2)
wg.to_html()

# %%
wg.run()


# %%
# Let's inspect the provenance graph for this workflow:

generate_node_graph(wg.pk)


# %%
# We can see that even though we passed the inputs as a single dictionary, they were serialized as two separate ``Int`` nodes, ``x`` and ``y``, before being passed to the task.
#
# Dynamic namespaces
# ~~~~~~~~~~~~~~~~~~
#
# Sometimes, the number and names of outputs are not known until the task runs.
# Dynamic namespaces are designed for this scenario.
#
# This is particularly useful for tasks that generate a variable number of outputs based on their inputs.

from aiida_workgraph import dynamic


@task
def generate_square_numbers(
    n: int,
) -> t.Annotated[
    dict,
    dynamic(t.Any),
]:
    """Generate a dict of square numbers. The number of outputs depends on 'n'."""
    return {f"square_{i}": i**2 for i in range(n)}


@task.graph
def generate_square_numbers_graph(n: int):
    generate_square_numbers(n=n)


wg = generate_square_numbers_graph.build_graph(n=5)
wg.run()

# %%
# Let's examine the provenance of this dynamic workflow:

generate_node_graph(wg.pk)

# %%
# The graph shows that the ``generate_square_numbers`` task has multiple output nodes, one for each entry in the dynamically generated dictionary.
# The ``dynamic(typing.Any)`` specification instructs the workgraph to treat each value in the returned dictionary as a separate output node of any type.
#
# Nested namespaces
# ~~~~~~~~~~~~~~~~~
#
# Namespaces can be nested to represent complex, structured data.
# Let's define a task that returns a nested dictionary.


@task
def nested_dict_task(
    x: int,
    y: int,
) -> t.Annotated[
    dict,
    namespace(sum=int, nested=namespace(diff=int, product=int)),
]:
    """Returns a nested dictionary with a corresponding nested namespace."""
    return {"sum": x + y, "nested": {"diff": x - y, "product": x * y}}


@task.graph
def nested_dict_graph(x: int, y: int):
    nested_dict_task(x=x, y=y)


wg = nested_dict_graph.build_graph(x=1, y=2)
wg.run()

# %%
# Instead of visualizing the full graph, let's inspect the outputs of the task using a summary utility.
# This is often clearer for verifying data structures.

print(get_process_summary(wg.tasks[-1].pk))

# %%
# The summary confirms that the output is correctly structured with a top-level ``sum`` and a nested ``nested`` dictionary, just as defined in the namespace.
#
# We can also combine dynamic and nested namespaces.


@task
def generate_dynamic_nested_dict(
    n: int,
) -> t.Annotated[
    dict,
    dynamic(namespace(square=int, cube=int)),
]:
    """Generate a nested dict of square and cube numbers from 0 to n."""
    return {f"data_{i}": {"square": i**2, "cube": i**3} for i in range(n)}


@task.graph
def generate_dynamic_nested_dict_graph(n: int):
    generate_dynamic_nested_dict(n=n)


wg = generate_dynamic_nested_dict_graph.build_graph(n=3)
wg.run()

# %%
# Let's check the output summary for this dynamically generated nested structure:

print(get_process_summary(wg.tasks[-1].pk))


# %%
# The output shows a dictionary with dynamic keys (``data_0``, ``data_1``, etc.), where each value is itself a dictionary with a fixed ``square`` and ``cube`` structure, as specified by ``dynamic(namespace(...))``.
#
# Data linkage
# ------------
#
# Data linkage tracks the flow of data between tasks.
# At the workflow level, a `task.graph` can define its own inputs and outputs, providing a clean interface to a complex chain of tasks.
# `aiida-workgraph` validates data against these graph-level specifications and automatically links graph inputs to the appropriate task inputs.
#
# In this final example, we will build a graph that reuses the input and output specifications from the tasks it contains.
# This is a powerful feature for building complex, modular, and self-consistent workflows.


@task
def add_multiply(
    data: t.Annotated[
        dict,
        namespace(x=int, y=int),
    ],
) -> t.Annotated[
    dict,
    namespace(sum=int, product=int),
]:
    """A reusable task with well-defined I/O specifications."""
    return {"sum": data["x"] + data["y"], "product": data["x"] * data["y"]}


@task.graph
def add_multiply_graph_final(
    n: int,
    data: t.Annotated[
        dict,
        namespace(
            add_multiply1=add_multiply.inputs,
            add_multiply2=add_multiply.inputs,
        ),
    ],
) -> t.Annotated[
    dict,
    namespace(
        square=generate_square_numbers.outputs,
        add_multiply1=add_multiply.outputs,
        add_multiply2=add_multiply.outputs,
    ),
]:
    """A complex graph demonstrating I/O reuse and data linkage."""
    square_numbers = generate_square_numbers(n)

    # Unpack nested inputs and pass them to the respective tasks
    out1 = add_multiply(data=data["add_multiply1"]["data"])
    out2 = add_multiply(data=data["add_multiply2"]["data"])

    # Gather task outputs into the graph-level output structure
    return {"square": square_numbers, "add_multiply1": out1, "add_multiply2": out2}


wg = add_multiply_graph_final.build_graph(
    n=3,
    data={
        "add_multiply1": {"data": {"x": 1, "y": 2}},
        "add_multiply2": {"data": {"x": 3, "y": 4}},
    },
)
# Generate an interactive HTML visualization of the graph
wg.to_html()

# %%
# In the example above:
#
# - **Graph outputs:** The ``outputs`` argument in the ``@task.graph`` decorator defines the *shape* of the final result.
#   We reuse ``generate_square_numbers.outputs`` and ``add_multiply_task.outputs`` to ensure the graph's output signature is consistent with the tasks it contains.
#
# - **Graph inputs:** The ``data`` input is annotated with a nested namespace that reuses ``add_multiply_task.inputs``.
#   This allows ``aiida-workgraph`` to validate the complex input dictionary and create the correct data links.
#
# In the GUI representation of the WorkGraph, you will see how the nested inputs are correctly wired.
# For instance, there is a direct link from the graph input socket ``data.add_multiply1.data.x`` to the task input socket ``add_multiply_task_1.data.x``, guaranteeing perfect data lineage.


# %%
# Conclusion
# ----------
#
# You now know how to annotate task and graph inputs and outputs in `aiida-workgraph`.
# By leveraging static (``namespace``), dynamic (``dynamic``), and nested namespaces, you can precisely control data serialization and create transparent data lineages.
#
# The key takeaways are:
#
# - Annotate task ``outputs`` to unpack results into individual AiiDA nodes.
# - Use ``typing.Annotated`` to specify input structures.
# - Employ ``dynamic`` for tasks with a variable number of outputs.
# - Reuse ``.inputs`` and ``.outputs`` specifications at the graph level to build modular and robust workflows.
#
# These tools are fundamental to building reproducible and verifiable scientific workflows with complete data provenance.
