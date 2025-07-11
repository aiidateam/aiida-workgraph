"""
=====================
How to use while loop
=====================
"""
# %%
# Introduction
# ============
# In this how-to we implement a while loop in ``WorkGraph``, allowing you to repeatedly execute a set of tasks as long as a specified condition remains true.
# We will implement the while logic using the while zone in the context manager and the modular task approach.
# To implement the while logic we need to use a context variable in both approaches, since we need a variable that can be reassigned in each iteration to compute the while-condition.
# This does not follow the regular `Dataflow programming <../../concept/autogen/dataflow_programming.py>`_ paradigm and therefore requires a specific construct like the context variable.
# If you are not familiar with context variables, please refer to `Use Context to pass data between tasks <context.ipynb>`_

# %%
# Setting up AiiDA environment
# ----------------------------

from aiida_workgraph import WorkGraph, task, While
from aiida_workgraph.utils import generate_node_graph
from aiida import load_profile

load_profile()


# %%
# While loop workflow
# ===================
# Suppose we want to implement a workflow with the following tasks:


@task
def compare(x, y):
    return x < y


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# The workflow logic expressed in Python looks like the following
#
# .. code-block:: python
#
#    n = add(1, 1)
#    while compare(n, 8):
#        n = add(n, 1)
#        n = multiply(n, 2)
#    result = add(n, 1)

# %%
# Context manager approach
# ------------------------
# Unlike regular tasks, the *while zone* lacks data input and output sockets.
# Tasks outside the zone can directly link to those inside, facilitating a workflow integration.
# We have the option to specify the maximum number of iterations to prevent an infinite loop.

with WorkGraph("while_context_manager") as wg:
    # set a context variable before running.
    outputs1 = add(x=1, y=1)
    wg.ctx.n = outputs1.result
    # We need to use context or static variables here
    outputs2 = compare(x=wg.ctx.n, y=8)
    # Because context variables do not follow the dataflow programming paradigm,
    # we have to state the dependency of the context variable explicitely.
    # The syntax below reads "`should_run` associated task waits for its executaion till `outputs1` has been successfully set its value"
    outputs2 << outputs1
    with While(outputs2.result, max_iterations=10):
        outputs3 = add(x=wg.ctx.n, y=1)
        outputs4 = multiply(x=outputs3.result, y=2)
        # We need to update the context variable since it is accessed in should_run task
        wg.ctx.n = outputs4.result
    wg.outputs.result = add(x=outputs4.result, y=1).result

wg.run()
print("State of WorkGraph:   {}".format(wg.state))
print("Result            :   {}".format(wg.outputs.result.value))
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view the ``While`` context manager is depicted as a *zone*, containing all its child tasks.
# This zone simplifies the visualization of the loop structure as it separates the logic executed within the loop from the one outside.
# Notice that the cyclic links around the context variable ``n`` since it is reused each iteration.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~

generate_node_graph(wg.pk)

# %%
# Task approach
# -------------
# Internally, the context manager approach uses the *while zone* task.
# We can also just directly add the *while zone* as a task using the ``workgraph.while_zone`` identifier.
# The task approach allows you to modularize your workflow separating the creation the workflow in different parts of your code.
# This allows you to create your code more dynamically, useful for an automatized creation of the workflow.

# ---

wg = WorkGraph("while_task")
initial_add_task = wg.add_task(add, x=1, y=1)
wg.ctx.n = initial_add_task.outputs.result

condition_task = wg.add_task(compare, x=wg.ctx.n, y=8)
condition_task.waiting_on.add(initial_add_task)
# --- while task starts
while_task = wg.add_task(
    "workgraph.while_zone", max_iterations=10, conditions=condition_task.outputs.result
)

add_task = while_task.add_task(add, x=wg.ctx.n, y=1)
multiply_task = while_task.add_task(multiply, x=add_task.outputs.result, y=2)
wg.ctx.n = multiply_task.outputs.result
# --- while task ends
wg.outputs.result = wg.add_task(add, x=multiply_task.outputs.result, y=1).outputs.result
wg.run()

print("State of WorkGraph:   {}".format(wg.state))
print("Result            :   {}".format(wg.outputs.result))
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15


# %%
# Workflow view
# ~~~~~~~~~~~~~
wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
generate_node_graph(wg.pk)


# %%
# Further reading
# ---------------
# Similarly, other the control logic can implemented, see `Use if condition  <../if.ipynb>`_ for if logic and `How to run tasks in parallel <parallel.py>`_ for map logic.
