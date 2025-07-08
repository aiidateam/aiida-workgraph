"""
==============================================
Graph Builder
==============================================
"""

# %%
# Introduction
# ============
# The `Graph Builder` is a powerful tool in the `aiida-workgraph` package that allows you to create dynamic workflows
# by defining a function that returns a `WorkGraph`. This function can be decorated with `@task.graph_builder` to
# specify the inputs and outputs of the workflow. The `Graph Builder` enables you to create workflows that can change
# depending on the inputs, making it particularly useful for integrating control flow logic such as loops and conditionals
# into your workflows.

# The `Graph Builder` is also the entry point for other people to use your workflows. User can call the function directly
# to create the `WorkGraph`.

# However, there are still some confusions about the `Graph Builder`
# - What is the use cases?
# - One can also use `If`, `Map` context manager to create dynamic workflows, why do we need the `Graph Builder`?
# - One can create nested workflows directly from a WorkGraph, why do we need the `Graph Builder`?
# - Another issue: how does the inputs of the `Graph Builder` serialize to AiiDA Data nodes? What data type should
#  I use for the inputs of the `Graph Builder`? AiiDA Data nodes or raw Python data types?

# In this document, we will answer these questions.

# Workflow factory
# ============================

# The `Graph Builder` is a function that is decorated with `@task.graph_builder`. The function returns a `WorkGraph`
# that defines the workflow. `Graph Builder` likes a factory function that creates a `WorkGraph` based on the inputs
# provided to the function.

from aiida_workgraph import task, WorkGraph
from aiida import orm, load_profile

load_profile()


@task()
def add(x, y):
    return x + y


@task.graph_builder()
def my_workflow(x, y):
    """ """
    wg = WorkGraph()
    wg.add_task(add, x=x, y=y)
    return wg


# Suppose you are a workflow developer and you created the your workflow, which is wrapped in the `my_workflow` function.
# For the users, they can use the `my_workflow` function to create a `WorkGraph` by passing the inputs directly as
# arguments. and then submit the workflow.

wg = my_workflow(1, 2)
wg.run()


# %%
# Control flow in workflows
# ============================

# The `Graph Builder` allows you to create dynamic workflows that can change depending on the inputs. This is particularly
# useful for integrating control flow logic such as loops and conditionals into your workflows.
# This method differs significantly from the `If` and `Map` context managers, which are used to create dynamic workflows.
#
# - **Visibility**: In the GUI, only the `graph_builder` task is visible before execution.
# - **Dynamic Generation**: Upon running, it generates the WorkGraph dynamically, allowing for complex conditional logic and flow adjustments based on runtime data.


@task()
def multiply(x, y):
    return x * y


@task.graph_builder(outputs=[{"name": "result"}])
def add_multiply_if(x, y):
    with WorkGraph() as wg:
        if x.value > 0:
            wg.outputs.result = add(x=x, y=y)
        else:
            wg.outputs.result = multiply(x=x, y=y)
        return wg


with WorkGraph() as wg:
    result = add(x=1, y=1)
    wg.add_task(add_multiply_if, x=result, y=2)
    wg.outputs.result = add(x=wg.tasks.add_multiply_if.outputs.result, y=1)


# %%
# Use `If` context manager
# -----------------------
# Use `If` context manager to create a dynamic workflow

from aiida_workgraph import task, WorkGraph, If

with WorkGraph() as wg:
    result = add(x=1, y=1)
    with If(result < 0):
        wg.ctx.data = add(x=result, y=2)
    with If(result >= 0):
        wg.ctx.data = multiply(x=result, y=2)
    # ---------------------------------------------------------------------
    result = add(x=wg.ctx.data, y=1)


# %%
# In the above example, it seems the `Graph Builder` need more lines of code to create the same workflow as the `If` context manager.
# However, inside the `Graph Buidler`, user can write normal Python code, instead of always need to use `Task`. Even through, the data provenance
# is not strictly keep. For example, if the first task return a dictionary, and the following task only use one of the value:


@task()
def sum_diff(x, y):
    return {"sum": x + y, "diff": x - y}


@task.graph_builder(outputs=[{"name": "result"}])
def add_multiply_if(x, y):
    with WorkGraph() as wg:
        if x.value["sum"] + 1 > 0:
            wg.outputs.result = add(x=x.value["sum"] + 1, y=y)
        else:
            wg.outputs.result = multiply(x=x.value["diff"] + 1, y=y)
        return wg


with WorkGraph() as wg:
    result = sum_diff(x=1, y=1)
    wg.add_task(add_multiply_if, x=result, y=2)
    wg.outputs.result = add(x=wg.tasks.add_multiply_if.outputs.result, y=1)

# If we use the context manager, we need to create another task to extract the value from the dictionary, which is not very convenient.


@task()
def extract_value(data, key):
    return data[key]


with WorkGraph() as wg:
    result = sum_diff(x=1, y=1)
    sum = extract_value(data=result, key="sum")
    condition = add(sum, 1)
    with If(condition < 0):
        data = add(sum, 1)
        wg.ctx.data = add(x=data, y=2)
    with If(condition >= 0):
        data = extract_value(data=result, key="diff")
        data = add(data, 1)
        wg.ctx.data = multiply(x=data, y=2)
    # ---------------------------------------------------------------------
    result = add(x=wg.ctx.data, y=1)


# So it's the balance between the convenience and the data provenance.
# The `Graph Builder` allows you to write more complex logic in a more convenient way,
# while the `If` context manager is more strict about the data provenance and the data flow
# in the workflow.
# Also Graph Builder allows users to write normal Python code, instead of domain specific language (DSL) like `If` context manager.


# %%
# Nested workflows with WorkGraph
# ===============================
# Both `Graph Builder` and `WorkGraph` allow you to create nested workflows. But what is the difference? and when should you use one over the other?

# The first use case is that you import the Graph Builder from another package, and you want to integrate it into your workflow.
# If all the inputs of the Graph Builder task are known at the time of the workflow creation, you can pass the inputs directly to the Graph Builder and create a `WorkGraph` directly.
# Then use the generated `WorkGraph` as a task in your workflow.
# If the inputs of the Graph Builder are not known at the time of the workflow creation, which means the input is the result of a previous task,
# then you need to use the `Graph Builder` to create the task.


#%%
# Data types for Graph Builder inputs
# ====================================

# As a developer, for user's convenience, the `Graph Builder` usually accepts raw Python data types as inputs, and the task will automatically serialize the data into AiiDA Data nodes during execution.
# However, the `Graph Builder` can also be used inside a `WorkGraph` to create nested workflows, in which case, the inputs of the `Graph Builder` task can be the results of previous tasks, which are AiiDA Data nodes.
# Because AiiDA track the provenance of the data, if the passed data is AiiDA data nodes, we should pass the data nodes to the tasks inside the `Graph Builder` directly, instead of using raw Python data types.
# Therefore, there is an inconsistency in the inputs of the `Graph Builder` task, which can be either raw Python data types or AiiDA Data nodes.
