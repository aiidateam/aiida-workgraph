"""
==================================================
Write workflows using the context manager paradigm
==================================================
"""

# %%
# This guide introduces the **context manager paradigm** in `aiida-workgraph`.
# The paradigm provides the means to explicitly define a workflow as a linear, non-nested sequence of tasks.
# This includes conditional and iterative constructs, allowing the user to clearly see the logical flow of the workflow.
# By using this approach to writing workflows, the user gains the power to control tasks not yet executed, as all possible paths are fully displayed.
#
# In the following sections, we'll explore how to build workflows using the context manager paradigm.
# We'll cover simple sequential workflows, conditional branches, and iterative execution in loops.
# Along the way, we'll highlight the contrasts with the ``@task``-based approach as they come up.
#
# Let's get started!

# %%
# Setup
# =====

from aiida import load_profile

load_profile()

# %%
# Creating a simple workflow
# ==========================
#
# We'll start by defining an add-multiply workflow.
# First, we define the necessary tasks:

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# Then the workflow itself. We first recall how this is done using the ``@task`` decorator:


@task.graph
def AddMultiply(x, y, z):
    the_sum = add(x, y).result
    return multiply(the_sum, z).result


# %%
# Now, let's see how we can achieve the same result using the context manager paradigm.

from aiida_workgraph import WorkGraph, spec
from typing import Any

with WorkGraph('AddMultiplyContextManager', inputs=spec.namespace(x=Any, y=Any, z=Any)) as wg:
    the_sum = add(
        x=wg.inputs.x,
        y=wg.inputs.y,
    ).result

    the_product = multiply(
        x=the_sum,
        y=wg.inputs.z,
    ).result

    wg.outputs.result = the_product


# %%
# A few things to note:
#
# - We explicitly name our workflow using the ``WorkGraph`` class
# - We get access to the ``WorkGraph`` instance (``wg``) directly
# - We use the ``WorkGraph`` context manager to define the workflow
# - Unlike the ``@task`` decorator, we explicitly set the inputs of the workflow using ``wg.inputs``
# - We pass to each task the inputs directly from the workflow's inputs
# - The results of the tasks are assigned to the workflow's outputs using ``wg.outputs.<socket_name>``
#
# How do they compare visually?

AddMultiply.build(x=1, y=2, z=3).to_html()

# %%

wg.to_html()

# %%
# With the context manager approach, we now gain a clear representation of workflow inputs and their connections to tasks.
# This is one benefit of explicitly defining workflow inputs and assigning input sockets to tasks explicitly.
#
# .. admonition:: Take home message
#
#    Using the context manager paradigm allows for clear visualization of workflow inputs and their connections to tasks.


# %%
# Workflow inputs/outputs
# =======================
#
# As we saw in the previous section, we get more control over the workflow inputs when using the context manager paradigm.
# This is provided via direct access to the ``WorkGraph`` instance.
# Let's see how we can use this further.

# %%
# Nested namespaces
# -----------------
#
# We can define our workflow inputs (and outputs) using namespaces for clarity and convenience.
# For example, consider the following workflow:

with WorkGraph(
    'AddThreeMultiplyContextManager',
    inputs=spec.namespace(
        add=spec.namespace(first=spec.namespace(x=Any, y=Any), second=spec.namespace(x=Any, y=Any)),
        multiply=spec.namespace(factor=Any),
    ),
    outputs=spec.namespace(sums=spec.namespace(first=Any, second=Any, third=Any), product=Any),
) as wg:
    first_sum = add(
        x=wg.inputs.add.first.x,
        y=wg.inputs.add.first.y,
    ).result

    second_sum = add(
        x=wg.inputs.add.second.x,
        y=wg.inputs.add.second.y,
    ).result

    third_sum = add(
        x=first_sum,
        y=second_sum,
    ).result

    product = multiply(
        x=third_sum,
        y=wg.inputs.multiply.factor,
    ).result

    wg.outputs = {
        'sums': {
            'first': first_sum,
            'second': second_sum,
            'third': third_sum,
        },
        'product': product,
    }

# %%
# When defining inputs under namespaces (here, ``add`` and ``multiply``), we can access them using the dot notation when assigning to tasks.
# Similarly, we can access the outputs using the same notation.
# Let's run our workflow, now with a clear input layout:


wg.run(
    inputs={
        'add': {
            'first': {
                'x': 1,
                'y': 2,
            },
            'second': {
                'x': 3,
                'y': 4,
            },
        },
        'multiply': {
            'factor': 5,
        },
    },
)

print('\nResults:')
print('  Sums:')
print('    First:', wg.outputs.sums.first.value)
print('    Second:', wg.outputs.sums.second.value)
print('    Third:', wg.outputs.sums.third.value)
print('  Product:', wg.outputs.product.value)

# %%
# Let's see what this looks like visually.

wg.to_html()

# %%
# Again, we see our workflow inputs, this time clearly organized under the ``add`` and ``multiply`` namespaces, with clear connections to each task.
# The same applies to the outputs, which are grouped under the ``sums`` and ``product`` namespaces.

# %%
# .. admonition:: Take home message
#
#    Input/output namespaces provide a clear and convenient way to organize and access workflow data.

# %%
# Input metadata
# --------------
#
# Workflow inputs can also be defined using ``wg.add_input(...)``.
# The method allows you to provide an identifier (e.g. ``workgraph.int``) to the input, which is used for validation.
#
# .. code:: python
#
#    with WorkGraph() as wg:
#        wg.add_input("workgraph.int", "x")  # validated as an integer
#        wg.add_input("workgraph.int", "y")
#        wg.add_input("workgraph.int", "z")
#        ...
#
# .. tip::
#
#    When using the AiiDA GUI, providing an ``identifier`` to input sockets will associate the input with a GUI component, allowing users to interact with the input in a more user-friendly type-specific way (see this :ref:`GUI section <web-ui:detailed-socket-view>`).
#
# In the future, further details may be added when defining inputs, e.g., default values, descriptions, help messages, etc.

# %%
# Nested workflows
# ================
#
# We can reuse existing workflows as tasks within our workflow.
# Let's have a look at how this works in practice:


@task
def generate_random_number(minimum, maximum):
    import random

    return random.randint(minimum, maximum)


def generate_add_multiply_workgraph():
    with WorkGraph(inputs=spec.namespace(x=Any, y=Any, z=Any), outputs=spec.namespace(result=Any)) as wg:
        the_sum = add(
            x=wg.inputs.x,
            y=wg.inputs.y,
        ).result

        the_product = multiply(
            x=the_sum,
            y=wg.inputs.z,
        ).result

        wg.outputs.result = the_product
    return wg


with WorkGraph(
    'AddMultiplyComposed',
    inputs=spec.namespace(min=Any, max=Any, x=Any, y=Any),
    outputs=spec.namespace(result=Any),
) as wg:
    random_number = generate_random_number(
        minimum=wg.inputs.min,
        maximum=wg.inputs.max,
    ).result

    nested_wg = generate_add_multiply_workgraph()

    wg.outputs.result = nested_wg(
        inputs={
            'x': wg.inputs.x,
            'y': wg.inputs.y,
            'z': random_number,
        }
    ).result

# %%
# We define a new workgraph, `AddMultiplyComposed` that reuses an `AddMultiply` workgraph (here wrapped in a reusable generator function) with a random number as input.
# In our new workflow, we first call the generator function to get an instance of the workgraph.
# We then call it with its inputs (similar to ``.run(inputs=...)``).
# As a task, it returns a socket namespace, in which we defined a ``result`` socket.
# Finally, we assign this ``result`` socket as the ``result`` socket of our composed workflow.
#
# Let's have a look at the graph:

wg.to_html()

# %%
# You can compare this graph against the one generated by the decorator paradigm (see :ref:`here <howto:combine-workgraphs:add-workgraph>`).
# Note the presence of the workflow inputs and their connections to the various tasks (as discussed earlier).

# %%
# Control flow
# ============
#
# In the decorator paradigm, we defer conditional and iterative logic to the ``@task.graph`` decorator to determine the flow at runtime.
# In the context manager paradigm, we instead define control flow logic explicitly with dedicated context managers.
# We explore this in the following sections.

# %%
# Workflow context
# ----------------
#
# Before we describe the control flow constructs, we need to understand the use of the workflow context.
# The details on AiiDA's context variables can be found in the :doc:`/advanced/autogen/context` advanced section.
# Here we describe how the context is treated by the ``WorkGraph``.
#
# ``WorkGraph`` doesn't track context variables (``wg.ctx``) for automatic dependency resolution because they could introduce cyclical dependencies between tasks.
# As such, when defining dynamic branching in a workflow (as is done in the following sections), we must explicitly do the following:
#
# 1. Store dynamically-generated results in the context variable
# 2. Define task dependencies using the ``<<`` or ``>>`` operators
#
# For further details, please refer to the :doc:`/howto/autogen/control_task_execution_order` how-to section.

# %%
# ``Zone``
# ---------
#
# A ``Zone`` acts as a single unit for dependency management, governed by two rules:
#
# 1. **Entry Condition**: A ``Zone`` (and all tasks within it) will only start
#    after *all* tasks with links pointing *into* the ``Zone`` are finished.
# 2. **Exit Condition**: Any task that needs an output from *any* task inside the
#    ``Zone`` must wait for the entire `Zone` to complete.
#
# This can be used to group a set of tasks.
# For example, consider the following workflow:


from aiida_workgraph import Zone


with WorkGraph('zone_example') as wg:
    task0_outputs = add(x=1, y=1)

    # This Zone will only start after task1 is finished,
    # because task3 depends on its result.
    with Zone() as zone1:
        task1_outputs = add(x=1, y=1)
        task2_outputs = add(x=1, y=task0_outputs.result)

    # Task 4 will wait for the entire Zone to finish,
    # even though it only needs the result from task2.
    task3_outputs = add(x=1, y=task1_outputs.result)

wg.to_html()

# %%
# The graph shows the first executed ``add`` task providing its result to the zone,
# then the two grouped ``add`` tasks inside the zone executing in parallel, and finally
# the last ``add`` task, which waits for the entire zone to finish before executing.

# %%
# ``If`` zone
# -----------
#
# To define the conditional logic of the workflow, ``WorkGraph`` provides an ``If`` context manager.
# Using the ``with If`` block, all child tasks are automatically encapsulated and executed only if the condition is met.

from aiida_workgraph import If
from aiida_workgraph.collection import group


with WorkGraph('AddMultiplyIf') as wg:
    first_sum = add(x=1, y=1).result

    with If(first_sum < 0):
        second_sum = add(x=first_sum, y=2).result
        wg.ctx.result = second_sum

    with If(first_sum >= 0):
        the_product = multiply(x=first_sum, y=3).result
        wg.ctx.result = the_product

    (final_sum := add(x=wg.ctx.result, y=1).result) << group(second_sum, the_product)

    wg.outputs.result = final_sum

wg.run()

print(f'Result: {wg.outputs.result.value}')
assert wg.outputs.result.value == 7

# %%
# Let's break it down:
#
# - We start by queuing an addition task
# - We then define the conditional branches using the ``If`` context manager
#
#   - We define the conditional logic w.r.t the ``result`` socket of the first task
#
#   .. note::
#
#      When sockets are involved in Pythonic expressions, they yield a task.
#      For example, ``first_sum < 0`` yields a built-in ``op_lt`` task.
#
#   - We define a task in each branch as we would do normally
#   - We assign the task result to the workflow context variable ``wg.ctx.result`` (see below for more details)
# - We define the dependency of ``final_sum`` on the result of the grouped,conditionally-determined tasks with the ``<<`` operator (``A << B`` means "A after B").

# %%
# .. note::
#
#    In the example above, as we are dealing with Python functions that are run in a blocking manner, the example would also work without explicit task dependency setting via ``<<`` or ``>>``, as further execution would anyway wait until both tasks have finished.
#    However, if the tasks would be submitted to the daemon in a non-blocking fashion (common use case in scientific scenarios with long-running jobs), the explicit waiting enforced by ``<<`` or ``>>`` is strictly required, so we also apply it here for consistency and correctness.
#
# Visualizing the graph, one can see two operator zones, ``op_lt`` and ``op_ge``, for our two comparisons ("less than" and "greater equal"), as well as one ``if_zone`` for each branch as defined by the two ``If`` context managers.
# Here, each ``if_zone`` has a ``conditions`` input socket, with both ``result`` sockets being fed into the ``graph_ctx``.
# From there, only one result is then fed as the input to the last add task (``add2``), and finally, the global workflow outputs.
# Lastly, we can see connections from each ``if_zone``'s special ``_wait`` output socket to the ``_wait`` input socket of the ``add2`` task, which represent the explicit waiting between the tasks as defined by the ``>>`` syntax.

wg.to_html()

# %%
# .. admonition:: Take home message
#
#    When using the ``If`` construct of the context manager paradigm, all branches are visible prior to execution.

# %%
# Finally, after the workflow has finished, we generate the provenance graph from the AiiDA process, where we can see that the result of the ``op_lt`` (less than) comparison is ``False`` and the branch ends there, while for the ``op_ge`` (greater or equal) comparison it is ``True``, meaning that the branch with the intermediate multiplication was executed.

wg.generate_provenance_graph()

# %%
# ``While`` zone
# --------------
#
# We can handle dynamically iterative tasks with the ``While`` context manager.
# Unlike regular tasks, the ``While`` zone lacks data input and output sockets.
# However, tasks outside the zone can directly link to those inside, facilitating workflow integration.
# We also have the option to specify the maximum number of iterations to prevent an infinite loop.

from aiida_workgraph import While


with WorkGraph('AddMultiplyWhile') as wg:
    n = add(x=1, y=1).result
    wg.ctx.n = n

    (condition := wg.ctx.n < 8) << n

    with While(condition, max_iterations=10):
        n = add(x=wg.ctx.n, y=1).result
        n = multiply(x=n, y=2).result
        wg.ctx.n = n

    wg.outputs.result = add(x=n, y=1).result

wg.run()

print(f'Result: {wg.outputs.result.value}')
# 2 -> While(3, 6 -> 7, 14) -> 15
assert wg.outputs.result.value == 15

# %%
# Once again, let's break it down:
#
# - We start by queuing an addition task to set the initial value of ``n``
# - We assign this initial value to the workflow context
# - We then define our condition, as well as the execution order (``condition`` waits for ``n``)
# - Next, we define the ``While`` context manager
#
#   - The block executes its inner tasks as long as the condition is met
#   - We define the inner tasks that will be executed iteratively (note the use of ``wg.ctx.n`` to access the value of ``n`` at the beginning of each iteration)
#   - We update the value of ``n`` in the workflow context at each iteration
# - Finally, we define the output of the workflow as the result of the last addition task
#
#   .. important::
#
#      Unlike the case of ``If``, the final addition task does not need to wait explicitly on the ``While`` zone, as the zone is guaranteed to be executed (unlike the two ``If`` branches).
#      Hence, the execution order is implicitly defined by the natural order of task definition.
#
# Visually, the ``While`` context manager is depicted as a `while zone`, containing all its child tasks.
# The zone simplifies the visualization of the loop structure, as it separates the logic executed within the loop from the one outside.

wg.to_html()

# %%
# .. note::
#
#    The cyclic links around the context variable ``n`` are due to its reuse in each iteration.
#
# In the provenance graph, we can see the looping and execution of multiple tasks in the loop reflected in the deep tree structure:

wg.generate_provenance_graph()

# %%
# ``Map`` zone
# ============
# .. warning::
#   **This feature is experimental.** The API for ``Map`` zone is subject to change in future releases.
#   We welcome your feedback on its functionality.
#   **ctx** does not work inside a ``Map`` zone yet.
#
# The ``Map`` context manager works similarly to Python's built-in ``map`` function.
# By accessing the ``item`` member of the ``Map`` context, we can pass each individual element (e.g. a dictionary entry) to tasks.
# This creates a new task behind the scenes for each element.

# %%
# Running tasks in parallel
# -------------------------
#
# One use case of ``Map`` is to run tasks in parallel.
# Let's see how this is done.
# We first define our data and a data extraction task:


len_list = 4
data = {f'data_{i}': {'x': i, 'y': i} for i in range(len_list)}


@task
def get_value(data, key):
    return data[key]


# %%
# .. note::
#
#    To perform an addition operation, we must extract the ``x`` and ``y`` values in separate tasks.
#    This is because tasks cannot be nested within other tasks - one of AiiDA's core tenets ensuring strict tracking of created data.
#
# Now we put it all together in a workflow that sums pairs of numbers in parallel:

from aiida_workgraph import Map


with WorkGraph('AddMap') as wg:
    with Map(data) as map_zone:
        result = add(
            x=get_value(map_zone.item.value, 'x').result,
            y=get_value(map_zone.item.value, 'y').result,
        ).result
        map_zone.gather({'result': result})
        wg.ctx.result = map_zone.outputs.result
        wg.outputs.result = wg.ctx.result

wg.run()

print('\nResults:')
for i, item in enumerate(wg.outputs.result.value.values()):
    print(f'  Item {i}: {item}')
# (1+1) + (2+2) + (3+3) = 12
assert sum(wg.outputs.result.value.values()) == 12

# %%
# Let's inspect the task graph:

wg.to_html()

# %%
# We can see with the ``Map`` zone the tasks onto which the input data would be mapped.
#
# .. tip::
#
#    If you are using the AiiDA GUI, you can visualize each instance of the mapped tasks by inspecting the ``Map`` zone.

# %%
# Finally, let's have a look at the provenance graph:

wg.generate_provenance_graph()

# %%
# Data aggregation
# ----------------
#
# Another use case of ``Map`` is to aggregate results.
# Commonly known as a gather, aggregate, or reduce operation, it is often used to automatically analyze or summarize the output of parallel computations.


@task
def aggregate_sum(data: spec.dynamic(Any)) -> int:
    return sum(data.values())


with WorkGraph('AddAggregate') as wg:
    with Map(data) as map_zone:
        added_numbers = add(
            x=get_value(map_zone.item.value, 'x').result,
            y=get_value(map_zone.item.value, 'y').result,
        ).result
        map_zone.gather({'result': added_numbers})
    wg.outputs.result = aggregate_sum(map_zone.outputs.result).result

wg.run()

print('\nResult:', wg.outputs.result.value)
assert wg.outputs.result.value == 12


# %%
# .. _advanced:context-manager:continue-workflow:
#
# Continuing a workflow
# =====================
#
# .. warning::
#   **This feature is experimental.** The API is subject to change in future releases. We welcome your feedback on its functionality.
#
# One of the features of ``WorkGraph`` is its ability to continue previous workflows.
# When a workgraph finishes its execution, it saves its state in the AiiDA process node.
# This allows you to rebuild the workgraph from the process and add new tasks to continue the workflow.

# %%
# Modifying inputs
# ----------------
#
# Let's run an add-multiply workflow with a hardcoded multiplication factor:


# In order to restart a workflow, the tasks should be importable from a module (i.e. not defined on-the-fly).
from aiida_workgraph.tasks.tests import add, multiply

with WorkGraph('AddMultiplyToBeContinued', inputs=spec.namespace(x=Any, y=Any)) as wg1:
    the_sum = add(
        x=wg1.inputs.x,
        y=wg1.inputs.y,
    ).result

    the_product = multiply(
        x=the_sum,
        y=3,
    ).result

    wg1.outputs = {
        'sum': the_sum,
        'product': the_product,
    }

wg1.run(
    inputs={
        'x': 1,
        'y': 2,
    },
)

print('\nResults:')
print(f'  Sum: {wg1.outputs.sum.value}')
print(f'  Product: {wg1.outputs.product.value}')

# %%
# Suppose we now want to rerun it with a different multiplication factor.
# Let's see how that's done:

with WorkGraph.load(wg1.pk) as wg2:
    wg2.name = 'AddMultiplyModified'
    wg2.restart()
    wg2.tasks.multiply.inputs.y = 4

wg2.run()

print('\nResults:')
print(f'  Sum: {wg2.outputs.sum.value}')
print(f'  Product: {wg2.outputs.product.value}')

# %%
# Note that the sum has not changed (the ``value``, but more importantly, the ``pk``, as it is the same node).
# The product, however, is the result of the calculation repeating with the new input, hence a brand new node.

# %%
# Adding new tasks
# ----------------
#
# Let's now pick up the previous workgraph and extend it by a second addition, leveraging the results of the previous work.
from node_graph.socket_spec import add_spec_field, SocketSpec

with WorkGraph.load(wg2.pk) as wg3:
    wg3.name = 'AddMultiplyContinued'
    wg3.add_input('workgraph.any', name='z')  # introduce a new input socket
    # also need to update the inputs spec
    wg3._inputs = add_spec_field(wg3._inputs, 'z', SocketSpec(identifier='workgraph.any'))
    wg3.restart()
    new_sum = add(
        x=wg3.tasks.multiply.outputs.result,
        y=wg3.inputs.z,
    ).result
    wg3.outputs.new_sum = new_sum

wg3.to_html()

# %%
print(f'State of WorkGraph : {wg3.state}')
print(f'State of add       : {wg3.tasks.add.state}')
print(f'State of multiply  : {wg3.tasks.multiply.state}')
print(f'State of new add   : {wg3.tasks.add1.state}')

# %%
# Note the ``PLANNED`` new addition task. Let's run it.
# Let's run it with the new input:

wg3.run(
    inputs={
        'z': 5,
    },
)

print('\nResults:')
print(f'  Sum: {wg3.outputs.sum.value}')
print(f'  Product: {wg3.outputs.product.value}')
print(f'  New sum: {wg3.outputs.new_sum.value}')

# %%
# Again, note that the previous data nodes are the same.
# Only the new addition task ran and created a new data node.
#
# The full provenance remains intact, including the original workflow:


wg3.generate_provenance_graph()

# %%
# Conclusion
# ==========
#
# In this section, using the context manager paradigm, you learned how to:
#
# - Use ``wg.inputs = {...}`` to define many inputs at once, or to define nested (namespaced) inputs to group related parameters
# - Use ``wg.add_input(...)`` to define a graph-level input and provide additional metadata (e.g., type validation)
# - Use graph-level inputs in tasks (``wg.inputs.<name>``)
# - Use ``wg.outputs.<name>`` to expose graph-level outputs from internal tasks
# - Use the ``If`` context manager to define conditional branches in a workflow
# - Use the ``While`` context manager to define iterative tasks in a workflow
# - Use the ``Map`` context manager to distribute data across instances of a set of tasks in parallel
# - Load an existing workgraph using ``WorkGraph.load(pk)``
# - Continue a workgraph by modifying inputs or adding new tasks
