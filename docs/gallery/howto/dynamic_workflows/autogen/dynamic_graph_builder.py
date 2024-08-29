"""
===============================================
Use Graph Builder to create a dynamic workflows
===============================================

"""

# %%
# Introduction
# ============
# The `Graph Builder` allow user to create also dynamic workflows that can change depending on the input.
#
# Load the AiiDA profile.


from aiida import load_profile

load_profile()


# %%
# Dynamic workflow dependent on the input
# =======================================
# Suppose to run a different task depending on the input. We want to run the add_one task if the number is below 2
# otherwise we want to run a modulo 2 task. For thet we require to create a context.

from aiida_workgraph import task, WorkGraph
from aiida.orm import Int


@task.calcfunction()
def add_one(x):
    return x + 1


@task.calcfunction()
def modulo_two(x):
    return x % 2


@task.graph_builder(outputs=[{"name": "result", "from": "context.out"}])
def add_modulo(i: Int):
    wg = WorkGraph()
    if i.value < 2:
        task = wg.add_task(add_one, x=i)
    else:
        task = wg.add_task(modulo_two, x=i)

    task.set_context({"result": "out"})
    return wg


wg = WorkGraph()
task1 = wg.add_task(add_modulo, i=Int(1))
task2 = wg.add_task(add_modulo, i=task1.outputs["result"])
wg.to_html()

# %%
# Running the workgraph.

wg.run()
print("Output of first task", task1.outputs["result"].value)  # 1 + 1 result
print("Output of second task", task2.outputs["result"].value)  # 2 % 2 result

# %%
# Plotting provenance

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)
