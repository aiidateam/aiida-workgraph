"""
==============================================
Aggregate dynamically sized and nested outputs
==============================================
"""
# %%
# Introduction
# ============
# In this tutorial, you will learn how to aggregate dynamically sized outputs by even by linking tasks or by using context.
# Then we discuss at the end how to deal with nested datatypes.
#
# Load the AiiDA profile.
#


from aiida import load_profile

load_profile()


# %%
# Using multi-linking for dynamic outputs
# =======================================
#


from aiida_workgraph import task, WorkGraph
from aiida.orm import Int


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


wg = WorkGraph("aggregate_by_multilinking")

aggregate_task = wg.add_task(aggregate, name="aggregate_task")

# we have to increase the link limit because by default workgraph only supports one link per input socket
# this is still an experimental feature that is why
aggregate_task.inputs["collected_values"].link_limit = 50

for i in range(2):  # this can be chosen as wanted
    generator_task = wg.add_task(generator, name=f"generator{i}", seed=Int(i))
    wg.add_link(
        generator_task.outputs["result"],
        aggregate_task.inputs["collected_values"],
    )

wg.to_html()


# %%
# Run the workgrhap

wg.submit(wait=True)
print("aggregate_task result", aggregate_task.outputs["result"].value)


# %%
# Plot provenance

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Mulitple dynamically sized outputs
# ----------------------------------
#

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

# we have to increase the link limit because by default workgraph only supports one link per input socket
# this is still an experimental feature that is why
aggregate_task.inputs["collected_ints"].link_limit = 50
aggregate_task.inputs["collected_floats"].link_limit = 50


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

wg.run()
print("aggregate_task int_sum", aggregate_task.outputs["int_sum"].value)
print("aggregate_task float_sum", aggregate_task.outputs["float_sum"].value)


# %%
# Plot provenance

generate_node_graph(wg.pk)

# %%
# Using context for dynamic outputs
# =================================


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


# For this use case it is more convenient to use the graph_builder as we can expose the context under a more convenient name.
@task.graph_builder(
    outputs=[{"name": "result", "from": "context.generated"}]
)  # this port is created by `set_context`
def generator_loop(nb_iterations: Int):
    wg = WorkGraph()
    for i in range(nb_iterations.value):  # this can be chosen as wanted
        generator_task = wg.add_task(generator, name=f"generator{i}", seed=Int(i))
        generator_task.set_context({"result": f"generated.seed{i}"})
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
# Run the workgrhap

wg.run()
print("aggregate_task result", aggregate_task.outputs["result"].value)


# %%
# Plot provenance

generate_node_graph(wg.pk)

# %%
# To support multiple dynamically sized outputs we can add another context and link it.


# %%
# Nested datatypes
# ================
# In principle, these methods can be also used for nested data types.
# However AiiDA does not support a nesting of orm types which happens for lists and dicts.
# It therefore tries to convert the nested data structures to native types when it becomes part of the provenance.
# For example an `orm.Dict` of `orm.Int`s is converted to an `orm.Dict` of built-in integers.


from aiida.orm import Dict

some_dict = Dict({"key": Int(5)})
print(some_dict.get_dict())
some_dict.store()  # store it in the database thereby it becomes part of the provenance
print(some_dict.get_dict())


# %%
# If it cannot convert the orm.Data type to a built-in type, an error message will be thrown


from aiida.orm import StructureData
import numpy as np

nonserializable_data = StructureData(cell=np.random.rand(3, 3))
some_dict = Dict({"key": nonserializable_data})

try:
    some_dict.store()
except Exception as err:
    print(err)

# %%
# One has to therefore remove the nestedness by using built-in types for the collections.
# In this example we use a dict instead of an orm.Dict


@task.calcfunction()
def compute_cell_volumes(strucs):
    return {"result": {key: struc.get_cell_volume() for key, struc in strucs.items()}}


dict_of_strucs = {
    f"struc{i}": StructureData(cell=np.random.rand(3, 3)) for i in range(3)
}
try:
    compute_cell_volumes(dict_of_strucs)
except Exception as err:
    print(
        err
    )  # Oops still the same error, but we now used a dict instead of an orm.Dict?

# %%
# The problem that AiiDA tries to automatic convert built-in types to orm.Data types.
# Even though we passed a dict, AiiDA automatically converts it to an orm.Dict.
# This feature allows us to not wrap everything to orm.Data types, but is here obstructive.
# To resolve this issue we have to specify the input arguments


@task.calcfunction()
def compute_cell_volumes(**kwargs):
    strucs = kwargs["strucs"]
    return {
        "result": {key: Float(struc.get_cell_volume()) for key, struc in strucs.items()}
    }


dict_of_strucs = {
    f"struc{i}": StructureData(cell=np.random.rand(3, 3)) for i in range(3)
}
cell_volumes = compute_cell_volumes(strucs=dict_of_strucs)

# %%
# Plotting the provenance graph for one of the cell volumes
generate_node_graph(cell_volumes["struc0"].pk)
