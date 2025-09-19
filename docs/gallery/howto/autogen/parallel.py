"""
Run tasks in parallel (Scatter-Gather)
=====================================
"""

# %%
# Introduction
# ------------
#
# A common pattern in scientific workflows is "scatter-gather," where a collection of inputs is processed in parallel (the "scatter" phase), and the results are then collected for a final processing step (the "gather" phase). This is a powerful way to parallelize work *within* a single, larger workflow.
#
# For example, applying a squaring operation `x²` to each number in a list `[x₁, x₂, ..., xₙ]`. Each squaring operation can be performed independently of the others.
#
# This how-to demonstrates how to implement the scatter-gather pattern to leverage this feature.

import typing as t
from aiida import load_profile
from aiida_workgraph import namespace, task, dynamic

load_profile()


# %%
# Scatter
# ---------------------
#
# The "scatter" phase involves creating and running multiple independent tasks. We can achieve this by simply iterating over our inputs within a `WorkGraph` and creating a task for each item.
#
# Let's define a simple task to square a number.


@task
def square(x: int) -> int:
    """Square an integer."""
    return x * x


# %%
# Define a helper task that generates a dictionary of numbers.
# This will serve as input to the scatter phase.


@task
def generate_numbers(
    n: int,
) -> t.Annotated[dict[str, int], namespace(data=dynamic(int))]:
    """Generate a dictionary of numbers from 1 to n."""
    return {'data': {f'number_{i + 1}': i + 1 for i in range(n)}}


# %%
# At first glance, one might try to write the scatter logic directly as follows:
#
# .. code:: python
#
#     data = generate_numbers(n=n).data
#     squares = {}
#     # Since these tasks have no dependencies on each other, they run in parallel.
#     for key, value in data.items():
#         squares[key] = square(x=value).result
#
# However, **this will not work as expected**.
#
# The reason is that ``data`` is not immediately available when constructing the graph.
# Instead, ``generate_numbers(n=n).data`` is a *future output*, a placeholder that will only be resolved at runtime.
#
# To correctly handle this, we must wrap the loop inside another task graph. This ensures that the graph engine knows how to schedule and parallelize the tasks.


@task.graph
def ParallelSquare(
    data: t.Annotated[dict[str, int], dynamic(int)],
) -> t.Annotated[dict, namespace(squares=dynamic(int))]:
    """Applies the square task to each number in the input data in parallel."""
    squares = {}
    for key, value in data.items():
        squares[key] = square(x=value).result
    return {'squares': squares}


# %%
# .. tip::
#
#     - **Parallelization**: In `WorkGraph`, tasks without data dependencies between them are automatically scheduled to run in parallel.
#     - **Dynamic inputs/outputs**: Our use of dynamic type annotations, such as `dynamic(int)`, allows AiiDA to create a distinct node for each input and output in the collection, which is essential for data provenance tracking. For more details, please refer to the section on :ref:`Dynamic namespaces <dynamic_namespaces>`.
#

# %%
# Let's run it with some sample data.

data = {f'number_{i}': i for i in range(1, 5)}

wg = ParallelSquare.build(data)
wg.run()

print('\nResults:')
for i, result_node in enumerate(wg.outputs.squares):
    original_value = list(data.values())[i]
    print(f'{original_value}² = {result_node.value}')

# %%
# Workflow view
# """""""""""""

wg.to_html()

# %%
# Provenance graph
# """"""""""""""""

wg.generate_provenance_graph()

# %%
# Gather
# --------------------
#
# The "gather" phase involves collecting the results from the parallel tasks and performing a final operation.
#
# We will now extend the workflow by adding a task that sums the results from the `square` tasks.


@task
def gather_and_sum(data: t.Annotated[dict, dynamic(int)]) -> int:
    """Sums the values of a dictionary of integers."""
    return sum(data.values())


# %%
# We create a new `WorkGraph` that orchestrates the full scatter-gather pattern.
# It first calls our `generate_numbers` and `ParallelSquare` graphs (scatter),
# and then feeds the collected outputs into the `gather_and_sum` task (gather).


@task.graph
def ScatterGatherSquare(n: int) -> int:
    """A full scatter-gather workflow to generate numbers, square them in parallel, and sum the results."""
    # Generate inputs
    data = generate_numbers(n=n).data
    # Scatter phase
    squares = ParallelSquare(data=data).squares
    # Gather phase
    return gather_and_sum(data=squares).result


wg = ScatterGatherSquare.build(4)
wg.run()

print('\nAggregated Result:', wg.outputs.result.value)

assert wg.outputs.result.value == 30

# %%
# Conclusion
# ----------
#
# In this how-to, we demonstrated how to implement the powerful scatter-gather pattern using `WorkGraph`.
#
# - The **scatter** phase is achieved by creating multiple independent tasks within a graph, which the engine automatically runs in parallel.
# - The **gather** phase collects the results from the parallel tasks for a final processing step.
# - We also highlighted an important concept: **future outputs**. When a value is the output of another task, it cannot be used in a Python loop directly at graph-construction time. Instead, the loop must be wrapped in another task graph.
#
