"""
=====================
Run tasks in parallel
=====================
"""

# %%
# Introduction
# ============
# Once you have developed a correct and functioning workflow, the next step is often to scale it up for large datasets.
# This typically involves applying the same workflow to many independent data points.
# In this how-to, we show how to run the workflow in parallel for each data point to improve performance and scalability.

# %%
# Setting up the AiiDA environment
# --------------------------------

from aiida import load_profile

load_profile()

from aiida_workgraph.utils import generate_node_graph
from aiida_workgraph import WorkGraph, task, spec

# %%
# Perfectly parallelizable problem
# ================================
# A perfectly parallelizable problem can be broken down into smaller, independent subproblems that require no shared resources.
# For example, consider an addition operation ``x + y`` applied element-wise to two lists: ``[x₁, ..., xₙ]`` and ``[y₁, ..., yₙ]``.
# Each individual addition can be performed independently of the others.
# ``WorkGraph`` automatically parallelizes task execution when there are no data dependencies between tasks (for more details on this concept, refer to `WorkGraph Engine <../../concept/autogen/engine>`_).
# We will take advanatge of this concept and create three different show three different ways how one can parallelize the add operation over the list with ``WorkGraph``.
#
# .. note::
#
#     In practice, a simple element-wise addition like this would typically be parallelized at a lower level, such as using NumPy vectorization or multithreading.
#     However, we use it here for illustrative purposes.
#     The concepts demonstrated in this guide can be applied to any workflow that is perfectly parallelizable.


# %%
# Conventional for-loop
# ---------------------


@task
def add(x, y):
    return x + y


len_list = 4
x_list = list(range(len_list))
y_list = list(range(len_list))
sums = []

wg = WorkGraph("parallel_for_loop")
for i in range(len_list):
    add_task = wg.add_task(add, x=x_list[i], y=y_list[i])
    sums.append(add_task.outputs.result)

wg.run()
print("Result:", sums)
# (1+1) + (2+2) + (3+3) = 12
assert sum(sum_socket.value for sum_socket in sums) == 12

# %%
# Workflow view
# ~~~~~~~~~~~~~
wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
generate_node_graph(wg.pk)

# %%
# Graph builder
# -------------
# We continue with creating the same task for the graph builder.
# We will perform the for loop within the graph builder task.
#
# .. note::
#
#     We pack the two lists into a dictionary to pass the data to a task, because ``aiida-core`` supports dynamically sized data structures only through dictionaries.
#     While lists are supported to some extent, their usage is limited to primitive types.

from typing import Any, Annotated


@task.graph(outputs=spec.namespace(result=spec.dynamic(Any)))
def parallel_add_workflow(data) -> dict:
    result = {}
    for i, item in enumerate(data.values()):
        outputs = add(x=item["x"], y=item["y"])
        result[f"sum_{i}"] = outputs.result
    return {"result": result}


len_list = 4
data = {f"list_{i}": {"x": i, "y": i} for i in range(len_list)}

wg = WorkGraph("parallel_graph_task")
wg.add_task(parallel_add_workflow, data=data)
wg.outputs.result = wg.tasks.parallel_add_workflow.outputs.result
wg.run()
print("Result:")
for socket in wg.outputs.result:
    print(socket._name, socket.value)
# (1+1) + (2+2) + (3+3) = 12
assert sum(wg.outputs.result._value.values()) == 12

# %%
# Workflow view
# ~~~~~~~~~~~~~
wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
generate_node_graph(wg.pk)


# %%
# Gather results
# ==============
# We now extend the workflow by adding a task that sums the intermediate results.
# This step is commonly known as a gather, aggregate, or reduce operation.
# It is often used to automatically analyze or summarize the output of the parallel computations.


# %%
# Graph builder
# -------------
# We will extend it the whole workflow only by the ``aggregate_sum`` task
@task
def aggregate_sum(data: Annotated[dict, spec.dynamic(Any)]):
    return sum(data.values())


@task.graph(outputs=spec.namespace(result=spec.dynamic(Any)))
def parallel_add_workflow(data) -> dict:
    result = {}
    for i, item in enumerate(data.values()):
        outputs = add(x=item["x"], y=item["y"])
        result[f"sum_{i}"] = outputs.result
    return {"result": result}


len_list = 4
data = {f"list_{i}": {"x": i, "y": i} for i in range(len_list)}

wg = WorkGraph("parallel_graph_task")
wg.add_task(parallel_add_workflow, data=data)
wg.add_task(aggregate_sum, data=wg.tasks.parallel_add_workflow.outputs.result)
wg.outputs.result = wg.tasks.aggregate_sum.outputs.result
wg.run()
print("Result:", wg.outputs.result.value)
assert wg.outputs.result == 12
