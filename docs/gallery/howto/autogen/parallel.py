"""
Run tasks in parallel
=====================
"""

# %%
# Introduction
# ------------
#
# Once you have developed a correct and functioning workflow, the next step is often to scale it up for large datasets.
# This typically involves applying the same workflow to many independent data points.
# In this how-to, we show how to run the workflow in parallel for each data point to improve performance and scalability.

import typing as t

from aiida import load_profile

from aiida_workgraph import namespace, task, dynamic


load_profile()

# %%
# Perfectly parallelizable problem
# --------------------------------
#
# A perfectly parallelizable problem can be broken down into smaller, independent subproblems that require no shared resources.
# For example, consider an addition operation ``x + y`` applied element-wise to two lists: ``[x₁, ..., xₙ]`` and ``[y₁, ..., yₙ]``.
# Each individual addition can be performed independently of the others.
# ``WorkGraph`` automatically parallelizes task execution when there are no data dependencies between tasks (for more details on this concept, refer to `WorkGraph Engine <../../concept/autogen/engine>`_).
#
# We will take advantage of this concept and show three different ways of how one can parallelize the add operation over the list with ``WorkGraph``.
#
# .. note::
#
#    In practice, a simple element-wise addition like this would typically be parallelized at a lower level, such as using NumPy vectorization or multithreading.
#    We use it here for illustrative purposes.
#    The concepts demonstrated in this guide can be applied to any workflow that is perfectly parallelizable.

# %%
# Conventional for-loop
# ~~~~~~~~~~~~~~~~~~~~~


@task
def add(x: int, y: int) -> int:
    return x + y


@task.graph
def ParallelAdd(
    data: t.Annotated[dict[str, dict[str, int]], dynamic(namespace(x=int, y=int))],
) -> t.Annotated[dict, namespace(sums=dynamic(int))]:
    sums = {}
    for i, list_i in enumerate(data.values()):
        x, y = list_i.values()
        sums[f"sum_{i}"] = add(x=x, y=y).result
    return {"sums": sums}


# %%
# Let's run it with some sample data.

data = {f"list_{i}": {"x": i, "y": i} for i in range(1, 5)}

wg = ParallelAdd.build(data)
wg.run()

print("\nResults:")
lists = list(data.values())
for i, sum_ in enumerate(wg.outputs.sums):
    print(f"{lists[i]['x']} + {lists[i]['y']} = {sum_.value}")

# %%
# Each addition was executed independently, yielding an AiiDA node for each result.

# %%
# Workflow view
# """""""""""""

wg.to_html()

# %%
# Provenance graph
# """"""""""""""""

wg.generate_provenance_graph()

# %%
# .. note::
#
#    Due to our explicit use of type annotation, AiiDA yields a node per input/output.
#    For more on leveraging type annotations, refer to the :doc:`/howto/autogen/annotate_inputs_outputs` how-to section.

# %%
# Gather results
# ~~~~~~~~~~~~~~
#
# We now extend the workflow by adding a task that sums the intermediate results.
# This step is commonly known as a gather, aggregate, or reduce operation.
# It is often used to automatically analyze or summarize the output of parallel computations.
#
# We will extend it the whole workflow only by the ``aggregate_sum`` task


@task
def aggregate_sums(data: t.Annotated[dict, dynamic(int)]) -> int:
    return sum(data.values())


@task.graph
def ParallelAddAggregate(
    data: t.Annotated[dict[str, dict[str, int]], dynamic(namespace(x=int, y=int))],
) -> int:
    sums = ParallelAdd(data=data).sums
    return aggregate_sums(data=sums).result


wg = ParallelAddAggregate.build(data)
wg.run()

print("\nResult:", wg.outputs.result.value)

# (1+1) + (2+2) + (3+3) + (4+4) = 20
assert wg.outputs.result.value == 20

# %%
# Conclusion
# ----------
#
# In this how-to, we demonstrated how to run tasks in parallel using ``WorkGraph``.
# We illustrated this with a simple example of element-wise addition, showcasing how to structure a perfectly parallelizable workflow.
# We also showed how to gather the results of parallel computations using an aggregation task.
