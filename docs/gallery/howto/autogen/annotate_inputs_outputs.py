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

from aiida_workgraph import dynamic, namespace, task
from aiida_workgraph.utils import get_process_summary

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


@task
def add_multiply1(x, y):
    """Return a dictionary, which will be stored as a single Dict node."""
    return {'sum': x + y, 'product': x * y}


@task
def add_multiply2(
    x: int,
    y: int,
) -> t.Annotated[
    dict,
    namespace(sum=int, product=int),
]:
    """Return a dictionary, but its elements are stored as separate Int nodes."""
    return {'sum': x + y, 'product': x * y}


@task.graph
def AddMultiply(x: int, y: int):
    """A graph to run both versions of the add_multiply task."""
    add_multiply1(x=x, y=y)
    add_multiply2(x=x, y=y)


wg = AddMultiply.build(x=1, y=2)
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

wg.generate_provenance_graph()

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
    return {'sum': data['x'] + data['y'], 'product': data['x'] * data['y']}


@task.graph
def AddMultiplyInputs(x: int, y: int):
    add_multiply3(data={'x': x, 'y': y})


wg = AddMultiplyInputs.build(x=1, y=2)
wg.to_html()

# %%
# Note how the ``x`` input is passed to ``data.x`` (and similarly for ``y``).
# This is due to the namespace specifications.

wg.run()


# %%
# Finally, we can inspect the provenance graph for this workflow:

wg.generate_provenance_graph()


# %%
# We can see that even though we passed the inputs as a single dictionary, they were serialized as two separate ``Int`` nodes, ``x`` and ``y``, before being passed to the task.
#
# .. _dynamic_namespaces:
# Dynamic namespaces
# ~~~~~~~~~~~~~~~~~~
#
# Sometimes, the number and names of outputs are not known until the task runs.
# Dynamic namespaces are designed for this scenario.
#
# This is particularly useful for tasks that generate a variable number of outputs based on their inputs.


@task
def generate_square_numbers(
    n: int,
) -> t.Annotated[
    dict,
    dynamic(t.Any),
]:
    """Generate a dict of square numbers. The number of outputs depends on 'n'."""
    return {f'square_{i}': i**2 for i in range(n)}


@task.graph
def SquareNumbersGenerator(n: int):
    generate_square_numbers(n=n)


wg = SquareNumbersGenerator.build(n=5)
wg.run()

# %%
# Let's examine the provenance of this dynamic workflow:

wg.generate_provenance_graph()

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
def generate_nested_dict(
    x: int,
    y: int,
) -> t.Annotated[
    dict,
    namespace(sum=int, nested=namespace(diff=int, product=int)),
]:
    """Returns a nested dictionary with a corresponding nested namespace."""
    return {'sum': x + y, 'nested': {'diff': x - y, 'product': x * y}}


@task.graph
def NestedDictGenerator(x: int, y: int):
    generate_nested_dict(x=x, y=y)


wg = NestedDictGenerator.build(x=1, y=2)
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
    return {f'data_{i}': {'square': i**2, 'cube': i**3} for i in range(n)}


@task.graph
def DynamicNestedDictGenerator(n: int):
    generate_dynamic_nested_dict(n=n)


wg = DynamicNestedDictGenerator.build(n=3)
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
# At the workflow level, a ``task.graph`` can define its own inputs and outputs, providing a clean interface to a complex chain of tasks.
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
    return {'sum': data['x'] + data['y'], 'product': data['x'] * data['y']}


@task.graph
def AddMultiplyFinal(
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
    out1 = add_multiply(data=data['add_multiply1']['data'])
    out2 = add_multiply(data=data['add_multiply2']['data'])

    # Gather task outputs into the graph-level output structure
    return {'square': square_numbers, 'add_multiply1': out1, 'add_multiply2': out2}


wg = AddMultiplyFinal.build(
    n=3,
    data={
        'add_multiply1': {'data': {'x': 1, 'y': 2}},
        'add_multiply2': {'data': {'x': 3, 'y': 4}},
    },
)
wg.to_html()

# %%
# In the example above:
#
# - **Graph outputs:** The outputs are annotated with a nested namespace that defines the *shape* of the final result.
#   Here we reuse ``generate_square_numbers.outputs`` and ``add_multiply_task.outputs`` to ensure the graph's output signature is consistent with the tasks it contains.
#
# - **Graph inputs:** The ``data`` input is annotated with a nested namespace that reuses ``add_multiply_task.inputs``.
#   This allows `aiida-workgraph` to validate the complex input dictionary and create the correct data links.
#
# In the GUI representation of the WorkGraph, you will see how the nested inputs are correctly wired.
# For instance, there is a direct link from the graph input socket ``data.add_multiply1.data.x`` to the task input socket ``add_multiply_task_1.data.x``, guaranteeing perfect data lineage.
#
# .. tip::
#
#    If a graph only exposes the outputs of a single task, this can be simplified as
#
#    .. code-block::
#
#       @task.graph
#       def SomeGraph(...) - t.Annotated[dict, some_task.outputs]:
#           return some_task(...)
#
# We can see similar linkage in the provenance graph.
# Let's run the graph and visualize its provenance.

wg.run()

wg.generate_provenance_graph()

# %%
# Note how the outputs of the various tasks are exposed (linked) to the graph, making accessible via the graph node.
#
# Reshaping specifications with ``select``
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Often, you'll want to reuse parts of a specification while modifying it.
# For example, you might want to exclude a field that is provided by another source or rename fields for clarity.
# This can be done declaratively using ``select``.
#
# The main parameters are:
#
# - ``include=...`` / ``exclude=...``: Keep or drop fields (supports **dotted paths** for nested fields).
# - ``include_prefix=...`` / ``exclude_prefix=...``: Filter top-level fields by name prefix.
# - ``rename={old:new}``: Rename top-level fields.
# - ``prefix="p_"``: Add a prefix to all top-level fields.
#
# .. Important::
#
#     To apply a ``select``, you must place it in the metadata list of ``t.Annotated[<type>, ...]`` alongside the specification you are modifying.
#
# Let's see an example where we build a graph that runs a task twice. We want a shared ``structure`` input at the graph level,
# so we must *exclude* it from the nested input specifications for each task call.

from aiida_workgraph.socket_spec import meta, select


@task
def consume_complex(
    data: t.Annotated[
        dict,
        namespace(
            pw=namespace(structure=int, kpoints=int, parameters=int),
            metadata=dict,
        ),
    ],
) -> dict:
    return {'seen': list(sorted(data.keys()))}


# Exclude a nested field (drop data.pw.structure)
# because it's provided separately as a graph input
# and shared between the two task calls.
@task.graph
def UseExclude(
    structure,
    inputs: t.Annotated[
        dict,
        namespace(
            consume_complex1=t.Annotated[dict, consume_complex.inputs, select(exclude='data.pw.structure')],
            consume_complex2=t.Annotated[dict, consume_complex.inputs, select(exclude='data.pw.structure')],
        ),
    ],
):
    # Manually add the shared 'structure' to each task's inputs
    consume_complex_input1 = inputs['consume_complex1']
    consume_complex_input1['data']['pw']['structure'] = structure
    consume_complex(**consume_complex_input1)

    consume_complex_input2 = inputs['consume_complex2']
    consume_complex_input2['data']['pw']['structure'] = structure
    consume_complex(**consume_complex_input2)


wg = UseExclude.build(
    structure=1,
    inputs={
        'consume_complex1': {'data': {'pw': {'kpoints': 2, 'parameters': 3}, 'metadata': {}}},
        'consume_complex2': {'data': {'pw': {'kpoints': 4, 'parameters': 5}, 'metadata': {}}},
    },
)
wg.to_html()

# %%
# In the GUI representation, you can see that the graph has a top-level ``structure`` input,
# and the nested inputs ``inputs.consume_complex1.data.pw`` and ``inputs.consume_complex2.data.pw`` are missing the ``structure`` socket, just as we specified.
#
# Modifying specification metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# You can also set namespace-level metadata declaratively using ``meta``, for example,
# to mark a reused specification as optional.
#
# - ``meta(required=False)``: Makes the input optional.
# - ``meta(is_metadata=True)``: Marks the input as metadata-only.
#
# It is attached alongside ``select`` in the same ``Annotated`` metadata list.


@task.graph
def UseMeta(
    data: t.Annotated[
        dict,
        consume_complex.inputs,  # Reuse the original spec
        meta(required=False),  # Make the entire 'data' input optional
    ],
):
    if data:
        return consume_complex(data=data)


# %%
# Conclusion
# ----------
#
# You now know how to annotate task and graph inputs and outputs in `aiida-workgraph`.
# By leveraging static (``namespace``), dynamic (``dynamic``), and nested namespaces, you can precisely control data serialization and create transparent data lineages.
#
# The key takeaways are:
#
# - Annotate task/graph outputs to unpack results into individual AiiDA nodes.
# - Annotate inputs to specify input structures.
# - Employ ``dynamic`` for tasks with a variable number of outputs.
# - Reuse ``.inputs`` and ``.outputs`` specifications at the graph level to build modular and robust workflows.
# - Use ``select`` inside ``Annotated`` to reshape reused specifications (e.g., ``include``/``exclude`` with dotted paths, ``rename``, `prefix`).
# - Use ``meta`` to modify metadata of a specification, such as making it optional.
#
# These tools are fundamental to building reproducible and verifiable scientific workflows with complete data provenance.
