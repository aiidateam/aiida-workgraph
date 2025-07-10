"""
===================================
Generate workflows programmatically
===================================
"""

# %%
# Introduction
# ============
# TODO

# %%
# Setting up AiiDA environment
# ----------------------------

from aiida_workgraph import WorkGraph, task, While
from aiida_workgraph.utils import generate_node_graph
from aiida import load_profile

load_profile()
# %%
# While loop
# ==========
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
# If logic
# ========
# Internally, the `if_` instruction is implemented using the `If Task` from the WorkGraph library. In the WorkGraph user interface, the `If Task` is visually represented as an "If Zone." This zone encapsulates all its child tasks, which are executed based on the defined conditions.
#
# - **Conditions**: The If Zone includes a `conditions` socket, which determines when the tasks inside the zone should be executed.
# - **Invert_condition**: If this input is True, it will invert the conditions.
# - **Task Linking**: Tasks located outside the If Zone can be directly linked to tasks within the zone, allowing for dynamic workflow adjustments based on conditional outcomes.
#
# Here is an example of how to add an `If Task` to a WorkGraph:
#
# ```python
# if2 = wg.add_task("If", name="if_false",
#                         conditions=condition1.outputs["result"],
#                         invert_condition=True)
# ```
#
# ### Adding tasks to the If Zone
# We can add tasks to the `If` zone using the `children` attribute.
#
# ```python
# # add task1 and task2 to the if zone
# if_task.children.add(["task1", "task2"])
# ```
#
# ### Creating the Workflow
# To construct the workflow, we'll utilize the built-in `If` and `Select` tasks from the Workgraph library. The `Select` task enables us to choose between two data sources based on a specified condition.
#
# The `Select task has the following inputs:
#
#    - **condition**: Provide the condition that dictates the selection between `true` and `false` outputs.
#    - **true**: Specify the output to be used if the condition evaluates to `true`. 
#    - **false**: Define the output for when the condition evaluates to `false`.
#

# %%
from aiida_workgraph import task, WorkGraph


with WorkGraph("if_task") as wg:
    condition = add(x=1, y=1)
    if_true_zone = wg.add_task("workgraph.if_zone", name="if_true",
                            conditions=condition)
    add2 = if_true_zone.add_task(add, name="add2", x=condition, y=2)
    if_false_zone = wg.add_task("workgraph.if_zone", name="if_false",
                            conditions=condition,
                            invert_condition=True)
    multiply1 = if_false_zone.add_task(multiply, name="multiply1", x=condition, y=2)
    #---------------------------------------------------------------------
    select1 = wg.add_task("workgraph.select", name="select1", true=add2.outputs["result"],
                        false=multiply1.outputs["result"],
                        condition = condition)
    add3 = wg.add_task(add, name="add3", x=select1.outputs["result"], y=1)
# export the workgraph to html file so that it can be visualized in a browser
wg.to_html()
# comment out the following line to visualize the workgraph in jupyter-notebook
# wg


# %%
# Map
# ===
# TODO
