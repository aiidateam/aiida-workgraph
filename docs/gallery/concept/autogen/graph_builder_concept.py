"""
==============================================
Graph Builder
==============================================
"""
# %%
# The *graph builder* provides a powerful and flexible way to create dynamic, reusable,
# and shareable workflows.
#
# At its core, a *graph builder* is a Python function decorated with
# ``@task.graph_builder`` that returns a ``WorkGraph`` instance. This design pattern allows
# you to use standard Python logic, including conditionals, loops, and complex data
# manipulation, to construct a workflow graph *dynamically* based on the inputs you provide.
#
# This guide addresses common questions about the *graph builder*:
#
# - What are its primary use cases?
# - How does it compare to context managers like ``If`` and ``Map``?
# - How does it differ from nesting a ``WorkGraph`` directly?
# - How should I handle input data types (Python vs. AiiDA nodes)?

# %%
# The `graph builder` as a workflow factory
# ============================================
#
# Think of a *graph builder* as a **factory for your workflows**. As a workflow
# developer, you encapsulate the logic for constructing a specific `WorkGraph` within a
# single, reusable function. Users can then call this function with their desired
# parameters to generate a ready-to-run `WorkGraph` instance without needing to
# understand its internal construction.
#
# This is the primary way you should expose your complex workflows to others.
#
# Example: A simple workflow factory
# ----------------------------------
#
# Here, `my_workflow` acts as a factory. A user can create and run a specific
# addition workflow by simply calling the function.

from aiida_workgraph import task, WorkGraph, If
from aiida import orm, load_profile

load_profile()


@task()
def add(x, y):
    return x + y


@task.graph_builder()
def my_workflow(x, y):
    """A simple workflow to add two numbers."""
    wg = WorkGraph("MyAdditionWorkflow")
    wg.add_task(add, x=x, y=y)
    return wg


# A user can now easily create and run the workflow:
wg = my_workflow(x=orm.Int(1), y=orm.Int(2))
# result = wg.run() # .run() is blocking, so we comment it out for docs.
print("Workflow created:", wg.name)


# %%
# Dynamic logic: `graph builder` vs. context managers
# ===================================================
#
# A key feature of the `graph builder` is enabling dynamic workflow structures. This
# often creates confusion when compared to using context managers like `If`, `While`,
# and `Map`. The fundamental difference lies in *when* and *how* the dynamic logic
# is executed.
#
# `graph builder`: pre-execution Python logic
# -------------------------------------------
#
# With a `graph builder`, the dynamic logic is **plain Python code that runs *before***
# **the AiiDA engine executes the workflow**. It constructs the `WorkGraph` based on
# the input values.
#
# **Advantages:**
#
# - **Convenience & Flexibility:** You can use the full power of Python to build the
#   graph. This is perfect for complex setup logic that doesn't need to be part of
#   the formal provenance graph.
#
# **Disadvantages:**
#
# - **Less Detailed Provenance:** The Python logic used to build the graph is not
#   tracked in the AiiDA provenance. This is a trade-off between convenience and
#   strict data provenance.
# - **Nesting Complexity:** A `graph builder` always creates a nested workflow. This can
#   complicate data flow, as inputs might need to be passed down through several
#   layers, whereas context managers operate within a "flat" graph.
# - **No Direct `while` Loop:** The `graph builder` pattern doesn't natively support
#   recurrent logic. For that, the `While` context manager is the appropriate tool.
#
# Example: conditional logic
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# In this example, the ``if data_node.value["sum"] + 1 > 0:`` is standard Python.
# It runs when the `add_multiply_if` graph builder is called, and only one branch
# (`add` or `multiply`) is ever added to the `WorkGraph`.


@task()
def multiply(x, y):
    return x * y


@task()
def sum_diff(x, y):
    return {"sum": x + y, "diff": x - y}


@task.graph_builder(outputs=[{"name": "result"}])
def add_multiply_if(data_node, y):
    """
    Builds a workflow that either adds or multiplies based on a
    value within the input data_node.
    """
    wg = WorkGraph()
    # Plain Python logic using the value of an AiiDA node
    if data_node.value["sum"] + 1 > 0:
        wg.outputs.result = add(x=data_node.value["sum"] + 1, y=y)
    else:
        wg.outputs.result = multiply(x=data_node.value["diff"] + 1, y=y)
    return wg


# --- Main workflow construction ---
with WorkGraph("GraphBuilderExample") as wg:
    sum_diff_result = sum_diff(x=1, y=1)
    # The add_multiply_if task will build and run its inner graph at execution time
    final_result = wg.add_task(add_multiply_if, data_node=sum_diff_result, y=2)


# %%
# `If` Context Manager
# ----------------------------------------------
#
# With the `If` context manager, the conditional logic is **part of the AiiDA**
# **`WorkGraph` itself**. Both branches of the condition exist in the graph, and the
# WorkGraph engine decides which path to execute based on the output of a preceding task.
#
# **Advantage:**
#
# - **Full Provenance:** Every logical step and data transformation is a node in the
#   graph, providing a complete and auditable record.
#
# **Disadvantages:**
#
# - **Verbosity:** Simple operations (like accessing a dictionary key) may require
#   dedicated `Task` nodes, making the workflow definition more verbose.
# - **Domain Specific Language (DSL):** You are constrained to the logic provided
#   by the context managers, not arbitrary Python.
#
# Example: the same logic with `If`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# To achieve the same result, we need an extra `extract_value` task. The workflow
# graph is flatter but more verbose.


@task()
def extract_value(data, key):
    return data[key]


# --- Main workflow construction ---
with WorkGraph("ContextManagerExample") as wg:
    result = sum_diff(x=1, y=1)
    sum_val = extract_value(data=result, key="sum")
    condition_val = add(sum_val, 1)

    with If(condition_val > 0):
        # This branch is defined in the graph
        data_to_use = add(sum_val, 1)
        wg.ctx.final_result = add(x=data_to_use, y=2)
    with If(condition_val <= 0):
        # This branch is also defined in the graph, using Else is idiomatic
        diff_val = extract_value(data=result, key="diff")
        data_to_use = add(diff_val, 1)
        wg.ctx.final_result = multiply(x=data_to_use, y=2)


# %%
# Key Differences at a Glance
# ---------------------------------
#
# +---------------------+--------------------------------------------+---------------------------------------------+
# | **Feature**         | **Graph Builder**                          | **Context Manager (`If`, `While`, `Map`)**  |
# +=====================+============================================+=============================================+
# | Flexibility         | High (any Python code for setup)           | Limited to the provided DSL                 |
# +---------------------+--------------------------------------------+---------------------------------------------+
# | Provenance          | Less detailed (setup logic is hidden)      | Fully detailed (logic is in the graph)      |
# +---------------------+--------------------------------------------+---------------------------------------------+
# | Structure           | Nested (creates a sub-workflow)            | Flat (adds nodes to the current graph)      |
# +---------------------+--------------------------------------------+---------------------------------------------+
# | Iterative Logic     | Not directly supported for ``while`` loops | Native support via ``While`` and ``Map``    |
# +---------------------+--------------------------------------------+---------------------------------------------+
# | Best For            | Reusable workflows, complex setup logic    | Strict provenance, iterative logic          |
# +---------------------+--------------------------------------------+---------------------------------------------+
#

# %%
# Nested workflows
# =====================
#
# Both a `graph builder` and a normal `WorkGraph` can be nested inside another
# `WorkGraph`. The choice depends on *when* the inputs to the sub-workflow are known.
#
# 1.  **Static Inputs (Known During Creation):**
#     If all inputs to your sub-workflow are known when you are building the main
#     workflow, you can call the `graph builder` function directly. This generates
#     a `WorkGraph` instance that you can then add as a single, nested task.
#

# Inputs (1, 2) are known upfront
sub_wg = my_workflow(1, 2)
main_wg = WorkGraph()
# Add the generated WorkGraph as a task
main_wg.add_task(sub_wg)


# %%
#
# 2.  **Dynamic Inputs (Result of a Previous Task):**
#     If the inputs to your sub-workflow depend on the output of another task in the
#     main workflow, you **must** use the `graph builder` as a task. The AiiDA engine
#     will wait for the inputs to be computed, then execute the `graph builder` to
#     generate and run the sub-workflow.
#

with WorkGraph() as main_wg:
    # 'add_result' is a future result from AiiDA.
    add_result = main_wg.add_task(add, x=10, y=5)
    # Use the Graph Builder as a task, feeding it the future result.
    main_wg.add_task(my_workflow, x=add_result, y=20)


# %%
# Handling Input Data Types
# ============================
#
# A `graph builder` must be robust enough to handle two scenarios for its inputs:
#
# 1.  **Raw python types:** When a user calls the `graph builder` function directly,
#     they will likely provide standard Python types (`int`, `str`, `dict`).
# 2.  **AiiDA data nodes:** When a `graph builder` is used as a task within a larger
#     `WorkGraph`, its inputs will be AiiDA data nodes (`orm.Int`, `orm.Str`,
#     `orm.Dict`) passed from previous tasks.
#
# Inside the graph builder's logic, you need to access the underlying Python value
# from an AiiDA node using its `.value` property to perform Python-native operations.
#
# Best Practice
# -------------
#
# Design your `graph builder` to handle both Python types and AiiDA nodes.
#
# **Note:** When passing the variable to tasks inside the `WorkGraph` (e.g.,
# `add(x=control_value, ...)`), you should pass the **original** variable
# (`control_value`), not the extracted Python value. This preserves the data
# provenance link if the input was an AiiDA node.


@task.graph_builder
def my_conditional_workflow(control_value, y):
    """
    A robust graph builder that handles Python types and AiiDA nodes.
    """
    wg = WorkGraph()

    # Use .value to safely access the data. If it is already a Python type,
    # this will raise an AttributeError, so we handle that case.
    try:
        py_value = control_value.value
    except AttributeError:
        py_value = control_value

    if py_value > 0:
        # Pass the original node/variable to preserve provenance
        wg.add_task(add, x=control_value, y=y)
    else:
        wg.add_task(multiply, x=control_value, y=y)

    return wg


# Example of using it within another workflow
with WorkGraph("RobustBuilderExample") as wg:
    result = add(x=-10, y=5)  # Result will be -5
    # The my_conditional_workflow will be executed with an AiiDA Int node
    wg.add_task(my_conditional_workflow, control_value=result, y=10)
