"""
==================================
Aggregate data from multiple tasks
==================================
"""
# %%
# Introduction
# ============
# In this tutorial, you will learn how to aggregate data from the outputs of
# multiple tasks by linking tasks or by using context. Then we discuss at the
# end how to deal with nested datatypes.
#
# By default, the workgraph only allows one link per input socket, mirroring
# typical function argument behavior where each argument accepts a single
# value. However, for scenarios requiring dynamic input port handling, e.g.,
# functions with variable keyword arguments, the link limit can be safely
# expanded to accommodate multiple inputs. An input port can be marked as dynamic
# by using a positional argument (`*args`) or a keyword argument (`**kwargs`) in the function. This adjustment is particularly
# useful in tasks designed to aggregate or process collections of data.


# Load the AiiDA profile.
from aiida import load_profile

load_profile()

# %%
# Dynamic input ports vs orm.Dict
# --------------------------
# We want to first focus shortly on the difference between a dynamic input ports and a orm.Dict as input
# as this is essential for understanding the aggregation mechanism. For that we firstly
# imitate the `aggregate` function using a Dict as input.

from aiida.orm import Int, Dict
from aiida_workgraph.utils import generate_node_graph
from aiida_workgraph import task, WorkGraph


@task.calcfunction()
def aggregate(
    **collected_values,
):  # We use the double asterisk to mark it as an dynamic input port
    for key, value in collected_values.items():
        print("key:", key, "value:", value)
    return Int(sum(collected_values.values()))


@task.calcfunction
def aggregate_dict(
    collected_values,
):  # We use the double asterisk to mark it as an dynamic input port
    for key, value in collected_values.items():
        print("key:", key, "value:", value)
    return Int(sum(collected_values.get_dict().values()))


some_dict = {f"value{i}": Int(i) for i in range(3)}
aggregate_sum = aggregate(**some_dict)
# Note that it is generally not recommended to have nested orm.Data types like
# in this case orm.Int in an orm.Data This will only work for json serializable
# types as orm.Int in this case. You will see in the provenance graph later
# that the information about the orm.Int's is completely lost.
aggregate_dict_sum = aggregate_dict(Dict(some_dict))

# %%
# We plot the provenance graph for both

generate_node_graph(aggregate_sum.pk)

# %%

generate_node_graph(aggregate_dict_sum.pk)

# %%
# We can see while the Dict version only considers the whole dictionary as one
# entity in the provenance graph, while the dynamic input port allows us to consider the different orm.Data instances.


# %%
# Using multi-linking to create dynamic input ports
# =================================================
# In the following example, we create multiple tasks that return a random
# integer, and then aggregate all the results and calculate the sum by linking
# the outputs of the tasks to the input of one final task


from aiida_workgraph import WorkGraph


@task.calcfunction()
def generator(seed: Int):
    import random

    random.seed(seed.value)
    return Int(random.randint(0, 1))  # one can also use


@task.calcfunction(outputs=[{"name": "sum"}])
def aggregate(
    **collected_values,
):  # We use the double asterisk to mark it as an dynamic input port
    for key, value in collected_values.items():
        print("key:", key, "value:", value)
    return {"sum": Int(sum(collected_values.values()))}


wg = WorkGraph("aggregate_by_multilinking")

aggregate_task = wg.add_task(aggregate, name="aggregate_task")

# we have to increase the link limit because by default workgraph only supports one link per input socket
aggregate_task.inputs["collected_values"].socket_link_limit = 50

for i in range(2):  # this can be chosen as wanted
    generator_task = wg.add_task(generator, name=f"generator{i}", seed=Int(i))
    wg.add_link(
        generator_task.outputs["result"],
        aggregate_task.inputs["collected_values"],
    )

wg.to_html()


# %%
# Run the workgraph

wg.submit(wait=True)

# %%
# Print the output
print("aggregate_task result", aggregate_task.outputs["sum"].value)


# %%
# The provenance is in this example also tracked. One can see how the generated
# integers are linked to the aggregate task and one final integer, the sum, is
# returned.

generate_node_graph(wg.pk)


# %%
# Multiple dynamic input ports
# ----------------------------
# We now do the same exercise as before but add another dynamic input port to it. We
# generate float numbers and link them to the aggregate task. The aggregate
# task now returns two sums one for the integers and one for the float numbers.
# To support this additional dynamic input ports, we can define multiple input
# sockets in the decorator, as shown in the code example below:


from aiida.orm import Float


@task.calcfunction()
def generator_int(seed: Int):
    import random

    random.seed(seed.value)
    return Int(random.randint(0, 1))


@task.calcfunction()
def generator_float(seed: Int):
    import random

    random.seed(seed.value)
    return Float(random.random())


# The variable keyword arguments (**collected_values) declare the input as dynamic.
# Thus, we can safely and flexibly add numerous input sockets.
@task.calcfunction(
    inputs=[{"name": "collected_ints"}, {"name": "collected_floats"}],
    outputs=[{"name": "int_sum"}, {"name": "float_sum"}],
)
def aggregate(**collected_values):
    print(collected_values)
    return {
        "int_sum": sum(collected_values["collected_ints"].values()),
        "float_sum": sum(collected_values["collected_floats"].values()),
    }


wg = WorkGraph("aggregate")

aggregate_task = wg.add_task(aggregate, name="aggregate_task")

# we have to increase the link limit because by default workgraph only supports
# one link per input socket.
aggregate_task.inputs["collected_ints"].socket_link_limit = 50
aggregate_task.inputs["collected_floats"].socket_link_limit = 50


for i in range(2):  # this can be chosen as wanted
    seed = Int(i)
    generator_int_task = wg.add_task(generator_int, name=f"generator_int{i}", seed=seed)
    generator_float_task = wg.add_task(
        generator_float, name=f"generator_float{i}", seed=seed
    )

    wg.add_link(
        generator_int_task.outputs["result"],
        aggregate_task.inputs["collected_ints"],
    )

    wg.add_link(
        generator_float_task.outputs["result"],
        aggregate_task.inputs["collected_floats"],
    )
wg.to_html()
# %%
# Run the workgraph

wg.submit(wait=True)

# %%
# Print the output
print("aggregate_task int_sum", aggregate_task.outputs["int_sum"].value)
print("aggregate_task float_sum", aggregate_task.outputs["float_sum"].value)


# %%
# Plot provenance

generate_node_graph(wg.pk)

# %%
# Using context for dynamic outputs
# =================================
# If your are not familiar with `context` please refer to the `doc page explaining it in detail <../context.html>`_.


@task.calcfunction()
def generator(seed: Int):
    import random

    random.seed(seed.value)
    return Int(random.randint(0, 1))  # one can also use


@task.calcfunction()
def aggregate(**collected_values):  # We use a keyword argument to obtain a dictionary
    for key, value in collected_values.items():
        print("key:", key, "value:", value)
    return {"result": sum(collected_values.values())}


# For this use case it is more convenient to use the graph_builder as we can
# expose the context under a more convenient name.
@task.graph_builder(
    outputs=[{"name": "result", "from": "context.generated"}]
)  # this port is created by `set_context`
def generator_loop(nb_iterations: Int):
    wg = WorkGraph()
    for i in range(nb_iterations.value):  # this can be chosen as wanted
        generator_task = wg.add_task(generator, name=f"generator{i}", seed=Int(i))
        generator_task.set_context({f"generated.seed{i}": "result"})
    return wg


wg = WorkGraph("generate_aggregate_by_context")


generator_loop_task = wg.add_task(
    generator_loop, name="generator_loop", nb_iterations=Int(2)
)

aggregate_task = wg.add_task(
    aggregate,
    name="aggregate_task",
    collected_values=generator_loop_task.outputs["result"],
)

wg.to_html()


# %%
# Run the workgraph

wg.submit(wait=True)

# %%
# Print the output

print("aggregate_task result", aggregate_task.outputs["result"].value)


# %%
# Plot provenance

generate_node_graph(wg.pk)

# %%
# To support multiple dynamically sized inputs we can add another context and link it.
