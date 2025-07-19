"""
=========================
Control flow in WorkGraph
=========================
"""
# %%
# Introduction
# ============
# In this how-to we show you how you can achieve the ``while`` loop and ``if`` conditional flow control elements in ``WorkGraph``.
# We'll explore two primary methods for achieving this:
#
# - The context manager approach (``If``, ``While``), which explicitly defines control flow zones within the graph.
# - The ``@task.graph`` decorator, which allows you to use native Python if/else statements and recursion to create dynamic, encapsulated workflows.
#
# So let's dive right into it!

# %%
# Setting up the AiiDA environment
# --------------------------------
# First, we need to set up our AiiDA environment, importing the necessary entities from ``aiida-workgraph`` and loading the AiiDA profile.

from aiida_workgraph import WorkGraph, task, If, While
from aiida_workgraph.utils import generate_node_graph
from aiida import load_profile

load_profile()


# %%
# If
# ==

# %%
# Workflow description
# --------------------
#
# Now, let's see how we can convert the following arithmetic workflow containing an ``if`` conditional into a WorkGraph:
#
# .. code-block:: python
#
#    # First addition
#    result = 1 + 1
#
#    # Conditionally execute addition or multiplication
#    if result < 0:
#        result = result + 2
#    else:
#        result = result * 3
#
#    # Last addition
#    result = result + 1

# %%
# First, we define the relevant arithmetic operations as WorkGraph tasks.
# Those will present the processes executed in the workflow, such that provenance is tracked.


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# Context manager
# ---------------
#
# To define the conditional logic of the workflow, WorkGraph provides an ``If`` context manager.
# Using the ``with If`` block, all the child tasks are automatically encapsulated.
# Their execution when the workflow is run is based on the defined conditions.

with WorkGraph("if_context_manager") as wg:

    add1_outputs = add(x=1, y=1)

    with If(add1_outputs.result < 0):
        cond_add_outputs = add(x=add1_outputs.result, y=2)
        wg.ctx.result = cond_add_outputs.result

    with If(add1_outputs.result >= 0):
        cond_mult_outputs = multiply(x=add1_outputs.result, y=3)
        wg.ctx.result = cond_mult_outputs.result

    final_outputs = add(x=wg.ctx.result, y=1)

    # -- Explicit dependency specification is required here --
    # Because wg.ctx.result is a context variable and not a direct task output,
    # WorkGraph cannot infer dependencies automatically (to avoid potential cycles).
    # Both branches (cond_add_outputs and cond_mult_outputs) may update wg.ctx.result,
    # so we must explicitly tell WorkGraph that ``final_outputs`` depends on each branch.
    # This is achieved with the ``<<`` and ``>>`` syntax.

    final_outputs << cond_add_outputs
    final_outputs << cond_mult_outputs

    # Finally, we set the result of the last task as the global workflow output result
    wg.outputs.result = final_outputs.result

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.outputs.result.value}")
assert wg.outputs.result.value == 7

# %%
# .. note:: WorkGraph doesn't track context variables (``wg.ctx``) for automatic dependency resolution because they could
#    introduce cyclical dependencies between tasks. By using the ``<<`` or ``>>`` operators, we explicitly declare that the
#    ``final_outputs`` task must wait for both ``cond_add_outputs`` and ``cond_mult_outputs`` to finish before running,
#    ensuring correct ordering when reading ``wg.ctx.result``. In the example above, as we are dealing with python
#    functions that are run in a blocking manner, the example would also work without explicit task dependency setting via
#    ``<<`` or ``>>``, as further execution would anyway wait until both tasks have finished. However, if the tasks would
#    be submitted to the daemon in a non-blocking fashion (common use case in scientific scenarios with long-running jobs),
#    the explicit waiting enforced by ``<<`` or ``>>`` is strictly required, so we also apply it here for consistency and
#    correctness.

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view, one can see two operator zones, ``op_lt`` and ``op_ge``, for our two comparisons
# ("less than" and "greater equal"), as well as one ``if_zone`` for each branch as defined by the two ``If`` context
# managers.
# Here, each ``if_zone`` has a ``conditions`` input socket, with both ``result``\ s being fed into the ``graph_ctx``.
# From there, only one result is then fed as the input to the last add task (``add2``), and, finally, the global ``graph_outputs``.
# Lastly, we can see connections from each ``if_zone``'s special ``_wait`` output socket to the ``_wait`` input socket of the ``add2`` task, which represent the explicit waiting between the tasks as request by the ``<<`` syntax.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# Finally, after the WG has finished, we generate the node (provenance) graph from the AiiDA process, where we can see
# that the result of the ``op_lt`` (larger than) comparison is ``False`` and the branch ends there, while for the ``og_ge`` (greater or equal) comparison it is ``True``, meaning that the branch with the intermediate multiplication was executed.

generate_node_graph(wg.pk)

# %%
# Graph Task
# ----------
#
# With the ``@task.graph`` decorator, we can create a Graph Task (see relevant section).
# This method differs significantly from the ``If`` context manager:
#
# - **Dynamic generation**: The WorkGraph is dynamically generated during runtime, allowing for complex conditional logic and flow adjustments based on runtime data.
# - **Visibility**: In the local workgraph view, only the ``graph`` task is visible before execution, with its internal
#   workings being hidden inside this *black box*. This is in contrast to the ``If`` context manager, for which both branches were shown.
# - **Use as task**: The *Graph Task* can be seamlessly added to other WGs, in the same way as a normal task, making the combination of multiple WGs easy.
# To achieve this, we use the ``@task.graph`` decorator, like so:

#%%
@task.graph()
def add_multiply_if(x, y, z):
    if x.value < 0:
        return add(x=x, y=y).result
    else:
        return multiply(x=x, y=z).result


#%%
# Inside the function body of our decorated function, we can thus write code in
# the same way as in a ``with WorkGraph`` context manager (that's one of the
# actual things the ``@task.graph`` decorators implicitly does).
# In this example, we define the three inputs, ``x``, ``y``, and ``z``, as we
# use different values for addition and multiplication, thus assigned to ``y``
# and ``z``.
# We can further directly use Python native ``if``, ``elif``, and ``else`` flow
# control elements.
# In each branch, we return the actual ``result`` of the respective task, which
# will be wrapped in the ``outputs`` (``TaskSocketNamespace``) of the *graph
# task*.

#%%
#
# Let's now see how we can use our graph task in another WG:

with WorkGraph("if_graph_task") as wg:
    first_add_result = add(x=1, y=1).result
    add_multiply_if_outputs = add_multiply_if(x=first_add_result, y=2, z=3)
    final_add_result = add(x=add_multiply_if_outputs.result, y=1).result

    wg.outputs.result = final_add_result

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.outputs.result.value}")
assert wg.outputs.result.value == 7

#%%
# As mentioned above, we can directly use our *graph task* the same as any
# other task, despite it containing an entire WG inside.
# Only its output namespace is returned, so we accessed the value via ``.result`` to pass it into the last task of the overall workflow.

# %%
# Workflow view
# ~~~~~~~~~~~~~
#
# In the local graphical workflow view, we only see the ``add_multiply_if1`` task, but the logic that was executed
# inside is hidden. Thus, the *graph task* presents somewhat of a "black box" in the graphical interface. This can be
# seen as a disadvantage, as it hides some of the actual workflow execution logic. Alternatively, it can also be seen as
# an advantage, for example in the case of complex top-level workflows that combine multiple sub-WorkGraphs, as it
# prevents a cluttering of the graphical interface by exposing all internal logic.  Instead in the WorkGraph web UI (see
# dedicated section), the task also appears as a black box  *before* execution, however, once it is executed, it can be
# expanded, and the individual tasks inside be visualized.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# In the provenance graph, we see that another WorkGraph was created and executed, as orchestrated by the top-level WG.
# Again, only the ``multiply`` branch was run with the given values.

generate_node_graph(wg.pk)

# %%
# While
# =====

# %%
# Workflow description
# --------------------
# Suppose we have an arithmetic workflow in Python with the following logic:
#
# .. code-block:: python
#
#    n = 1 + 1
#    while n < 8:
#        n = n + 1
#        n = n * 2
#        print(n)
#    result = n + 1

# %%
# To convert this simple workflow into a WorkGraph, we again require the necessary ``task``\ s.
# As we already have the ``add`` and ``multiply`` tasks defined above, we only require one for the comparison:


@task
def compare_lt(x, y):
    return x < y


# %%
# Context manager
# ---------------
# To construct the WorkGraph, we can again use the ``While`` context manager,
# which effectively creates a ``while zone`` (as seen in the workflow view below).
# Unlike regular tasks, the *while zone* lacks data input and output sockets.
# But, tasks outside the zone can directly link to those inside, facilitating workflow integration.
# Finally, we have the option to specify the maximum number of iterations to prevent an infinite loop.
# The whole code snippet is shown here, and additional comments are given below:

with WorkGraph("while_context_manager") as wg:
    initial_n = add(x=1, y=1).result
    wg.ctx.n = initial_n
    should_run = compare_lt(x=wg.ctx.n, y=8).result
    should_run << initial_n
    with While(should_run, max_iterations=10):
        n = add(x=wg.ctx.n, y=1).result
        n = multiply(x=n, y=2).result
        wg.ctx.n = n

    wg.outputs.result = add(x=n, y=1).result

wg.run()
print(f"State of WorkGraph: {wg.state}")
print(f"Result            : {wg.outputs.result.value}")
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15

#%%
# First, we set ``initial_n`` as a context variable ``n`` before running the loop.
# We further define condition ``should_run`` using the ``compare_lt`` task.
# The syntax reads "``should_run`` associated task waits for its execution till ``outputs1`` has been successfully set its value" (see note above).
# These preparatory tasks set the stage (that is, create the necessary tasks, sockets, and links) in the WG, that we can
# now introduce the while-loop.
# This is achieved with the ``While`` context manager, in which the context variable ``n`` is continuously updated by
# the ``add`` and ``multiply`` operations.
# Due to the previously created links, this is reflected in the ``should_run`` task.
# Lastly, we execute the final addition once the while loop concludes.

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the graphical workflow view the ``While`` context manager is depicted as a *while zone*, containing all its child tasks.
# This zone simplifies the visualization of the loop structure as it separates the logic executed within the loop from the one outside.
# Notice that the cyclic links around the context variable ``n`` since it is reused each iteration.

wg.to_html()

# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# In the provenance graph, we can see the looping and execution of multiple tasks in the loop reflected in the deep tree structure:

generate_node_graph(wg.pk)

# %%
# Graph Task
# ----------
# We can also implement the same while loop logic using a graph task with the @task.graph decorator.
# Instead of a native while loop, this approach uses recursion: the graph task calls itself repeatedly until a termination condition is met.
# Each recursive call dynamically generates a new sub-workflow, effectively creating one "iteration" of the loop.
# First, let's define the recursive graph task.


@task.graph()
def add_multiply_while(n, N):
    """A recursive graph task that mimics a while loop."""
    # When n >= N, the recursion stops, and the current value of n is returned.
    if n >= N:
        return n
    n = add(x=n, y=1).result
    n = multiply(x=n, y=2).result
    # Call the function itself with the updated value of n.
    # This continues the loop, creating a new nested workflow.
    return add_multiply_while(n=n, N=N).result


# Now, we can use this recursive graph task within our main WorkGraph.
with WorkGraph("while_graph_task") as wg:
    first_add_result = add(x=1, y=1).result
    add_multiply_while_result = add_multiply_while(n=first_add_result, N=8).result
    final_add_result = add(x=add_multiply_while_result, y=1).result
    wg.outputs.result = final_add_result
    wg.run()
    print(f"Result            : {wg.outputs.result.value}")
    assert wg.outputs.result.value == 15


#%%
# This approach encapsulates the entire loop within a single, reusable task.

# %%
# Workflow view
# ~~~~~~~~~~~~~
# In the workflow view, add_multiply_while appears as a single "black box" task.

wg.to_html()


# %%
# Provenance graph
# ~~~~~~~~~~~~~~~~
# The provenance graph clearly show the chain of nested WorkGraph calculations, revealing each "iteration" as a distinct sub-process spawned by the recursive calls.

generate_node_graph(wg.pk)

# %%
# Conclusion
# ==========
#
# This tutorial demonstrates how to implement control flow structures (``if`` conditionals and ``while`` loops)
# in WorkGraph. The key concepts covered:
#
# - **If conditionals** can be implemented using either:
#
#   - The ``If`` context manager for explicit workflow visualization with visible branches
#   - The ``@task.graph`` decorator for dynamic runtime generation with encapsulation
#
# - **While loops** use:
#
#   - The ``While`` context manager to create iterative workflows with configurable maximum iterations to prevent infinite loops
#   - The ``@task.graph`` decorator, where loops are created through recursion to handle dynamic iterations.
#
# - **Context variables** (``wg.ctx``) require explicit dependency management using ``<<`` and ``>>`` operators
#   since WorkGraph cannot automatically infer dependencies to avoid potential cycles
# - **Provenance tracking** is maintained throughout all control flow operations during workflow execution
