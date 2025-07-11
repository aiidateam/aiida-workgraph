"""
============================
How to run tasks in parallel
============================
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
from aiida_workgraph import WorkGraph, task, Map

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


@task.graph_builder
def parallel_add_workflow(data):
    wg = WorkGraph()
    for i, item in enumerate(data.values()):
        add_task = wg.add_task(add, x=item["x"], y=item["y"])
        wg.outputs.result = {f"sum_{i}": add_task.outputs.result}
    return wg


len_list = 4
data = {f"list_{i}": {"x": i, "y": i} for i in range(len_list)}

wg = WorkGraph("parallel_graph_builder")
wg.add_task(parallel_add_workflow, data=data)
wg.outputs.result = wg.tasks.parallel_add_workflow.outputs.result
wg.run()
print("Result:", wg.outputs.result.value)
# (1+1) + (2+2) + (3+3) = 12
assert sum(wg.outputs.result.value.values()) == 12

# %%
# Workflow view
# ~~~~~~~~~~~~~
wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
generate_node_graph(wg.pk)

# %%
# Map context manager
# -----------
# The ``Map`` works similar as python's inbuilt map.
# By accessing the member ``item`` of the map context we can directly pass the socket item to tasks passing creating for each element a new task behind the curtain.
# There is a caveat, to apply the add operation we need to access the ``x`` and ``y`` elements in a separate task since we cannot run a task within a task.
# The ``Map`` context works similarly to Python's built-in ``map``.
# By accessing the ``item`` member of the ``Map`` context, we can pass each individual element (e.g. a dictionary entry) to tasks.
# This creates a new task behind the scenes for each element.
#
# .. note::
#
#   To perform an addition operation, we must extract the `x` and `y` values in separate tasks.
#   This is because tasks cannot be nested within other tasks, one of ``aiida-core`` concepts to be able to strictly track created data.


@task
def get_value(data, key):
    return data[key]


len_list = 4
data = {f"data_{i}": {"x": i, "y": i} for i in range(len_list)}

with WorkGraph("parallel_map") as wg:
    with Map(data) as map_:
        wg.outputs.result = add(
            x=get_value(map_.item, "x").result, y=get_value(map_.item, "y").result
        ).result

wg.run()
print("Result:", wg.outputs.result.value)
# (1+1) + (2+2) + (3+3) = 12
assert sum(wg.outputs.result.value.values()) == 12

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
def aggregate_sum(data):
    return sum(data.values())


@task.graph_builder
def parallel_add_workflow(data):
    wg = WorkGraph()
    for i, item in enumerate(data.values()):
        add_task = wg.add_task(add, x=item["x"], y=item["y"])
        wg.outputs.result = {f"sum_{i}": add_task.outputs.result}
    return wg


len_list = 4
data = {f"list_{i}": {"x": i, "y": i} for i in range(len_list)}

wg = WorkGraph("parallel_graph_builder")
wg.add_task(parallel_add_workflow, data=data)
wg.add_task(aggregate_sum, data=wg.tasks.parallel_add_workflow.outputs.result)
wg.outputs.result = wg.tasks.aggregate_sum.outputs.result
wg.run()
print("Result:", wg.outputs.result.value)
assert wg.outputs.result == 12

# %%
# Map context
# -----------
# Similarly for the map context approach we only need do extend it by the ``aggregate_sum`` task


@task
def aggregate_sum(data):
    return sum(data.values())


@task
def get_value(data, key):
    return data[key]


len_list = 4
data = {f"data_{i}": {"x": i, "y": i} for i in range(len_list)}

with WorkGraph("parallel_map") as wg:
    with Map(data) as map_:
        added_numbers = add(
            x=get_value(map_.item, "x").result, y=get_value(map_.item, "y").result
        ).result
    wg.outputs.result = aggregate_sum(added_numbers).result

wg.run()
print("Result:", wg.outputs.result.value)
assert wg.outputs.result == 12
