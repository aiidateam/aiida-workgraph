"""
============================
Wait for another task
============================
"""

# %%
# ## Introduction
# In WorkGraph, tasks are generally launched when their input tasks are completed.
# However, there may be scenarios where you need to delay the execution of a task until another unrelated task is complete.
# This can occur, for example, when the task's input depends on data dynamically populated in the `context`, and you must ensure this data is ready before proceeding.
#
# This tutorial will guide you through using the `>>` and `<<` operators to manage such dependencies effectively.
#
# .. note::
#
#    For cross-WorkGraph dependencies, i.e. you want to wait for a task from another WorkGraph, you can use the `monitor` task, please refer to the `Monitor <./monitor.rst>`__ tutorial.
#
#
#
# ## Example
# Here we create two `add` tasks and one `sum` task. The `sum` task will wait for the two `add` tasks to finish.

# %%
from aiida_workgraph import task
from aiida.orm import Float
from aiida import load_profile

load_profile()

# define add task
@task()
def add(x, y):
    return x + y


# define sum task
@task()
def sum(**datas):
    total = 0
    for data in datas.values():
        total += data.value
    return Float(total)


# %%


@task.graph()
def test_wait(x, y):
    from aiida_workgraph import get_current_graph

    wg = get_current_graph()
    wg.ctx.data = {}  # initialize context data
    outputs1 = add(x=x, y=1)
    wg.ctx.data.add1 = outputs1.result
    outputs2 = add(x=y, y=2)
    wg.ctx.data.add2 = outputs2.result
    # let sum task wait for add1 and add2, and so the `data` in the ctx is ready
    outputs3 = sum(datas=wg.ctx.data)
    outputs1 >> outputs3
    outputs2 >> outputs3
    return outputs3.result


# %%
# ### Run and check results
#
wg = test_wait.build_graph(x=1, y=2)
wg.run()
# %%
print("State of WorkGraph         : {}".format(wg.state))
print("Result of sum1: {}".format(wg.outputs.result.value))

# %%
# Generate node graph from the AiiDA process:

# %%
from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)
