"""
Quick Start
===========

"""


# %%
# Introduction
# ------------
#
# In this quick-start tutorial, we will cover:
#
# - the basics of ``WorkGraph`` with a simple arithmetic workflow
# - how to define and use
#   - tasks to encapsulate a process
#   - sub-workflows as tasks to encapsulate several sub-processes
# - an advanced example highlighting briefly the various features of ``WorkGraph``
# - a brief introduction to the web interface

# %%
# Installation
# ------------
#
# Let's first install ``aiida-workgraph``

# ::
#    pip install aiida-workgraph
#
# .. note::
#
#    This will also install ``aiida-core`` and its dependencies.
#
# .. note::
#    If this is your first time using AiiDA, you will need to configure a profile.
#    You can do this by running the following command in your terminal:

#    ::
#       verdi presto

# %%
# Load an AiiDA profile
# ---------------------
#
# To interact with the AiiDA database, you need to load your AiiDA profile.
# You have two options:
#
# 1. AiiDA terminal - You can run the following command to enter an AiiDA-configured ``IPython`` shell
#
#    ::
#       verdi shell
#
#    .. note::
#
#       This will start an IPython shell with your AiiDA profile loaded.
#       This environment is configured to work with AiiDA seamlessly, including such features as tab completion.
#
# 2. Alternatively, if working in a Jupyter notebook, you can load your profile with the following code:

from aiida import load_profile

load_profile()

# %%
# .. note::
#
#    You can find out more about AiiDA and its architecture in the `AiiDA documentation`_.

# %%
# (Optional) Start the AiiDA GUI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Some features of ``WorkGraph`` are best demonstrated interactively.
# For this, you can use the AiiDA GUI.
# First, let's install the necessary packages:
#
# ::
#
#    pip install aiida-gui aiida-gui-workgraph
#
# Once installed, open a terminal, and run:
#
# ::
#
#    aiida-gui start
#
# Go to http://localhost:8000/workgraph in your browser to view a table of executed ``WorkGraph`` processes.

# %%
# Simple workflow
# ---------------
#
# Suppose we want to calculate ``(x + y) * z`` in two steps:
#
# - add ``x`` and ``y``
# - then multiply the result with ``z``.
#
# The simplest way to do this in ``WorkGraph`` is as follows:

from aiida_workgraph import WorkGraph

with WorkGraph("AddMultiply") as wg:
    wg.inputs = dict.fromkeys(["x", "y", "z"])
    wg.outputs.result = (wg.inputs.x + wg.inputs.y) * wg.inputs.z

# %%
# We've defined our first workflow using ``WorkGraph`` 🎉
# It takes three arguments, ``x``, ``y``, and ``z``, adds the first two, multiplies the result by the third, and assigns the result as an output.
# In the background, ``WorkGraph`` converts the addition and subsequent multiplication operations into tasks.
#
# ..tip::
#
#   ``result`` is a special output port given by default.
#    However, you can also define your own output ports with ``wg.outputs.my_output = ...``.
#    You can learn more about graph-level inputs and output in the :doc:`autogen/graph_level` section.
#
# We can visualize the workgraph to inspect the tasks and their connections.

wg.to_html()

# We can run it with

wg.run(inputs={"x": 2, "y": 3, "z": 4})

# %%
# .. tip::
#
#    If you've configure a message broker in AiiDA, you can also submit the workgraph
#    (queue it for execution) with ``wg.submit(inputs=...)``.
#    You can learn more about this in the `AiiDA documentation_`.

# %%
# .. tip::
#
#    If you're running the AiiDA GUI, you can visualize the executed workflow interactively.
#    Click on the ``PK`` field of the submitted workflow (look for ``WorkGraph<AddMultiply>``) to view its details.

# %%
print("Result:", wg.outputs.result.value)

# %%
# We can also access the results of the individual tasks.

print("Result of addition:", wg.tasks.op_add.outputs.result.value)
print("Result of multiplication:", wg.tasks.op_mul.outputs.result.value)

# %%
# .. note::
#
#    Many Python operators are supported in `WorkGraph`` and are automatically converted to tasks with the name ``op_<operator>``, e.g. ``op_add``, ``op_mul``, ``op_lt``, ``op_eq``, and more.
#
# So far so good! We have a simple workflow, can visualize it, run/submit it with inputs, and inspect its outputs.
# However, in practice, tasks will require operations more complex than simple arithmetics.
# How can we handle these?

# %%
# ``Task``s
# ---------
#
# A ``Task`` is the basic building block of the ``WorkGraph``. It has inputs,
# outputs, and an executor. Let's define the above operations as tasks.

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# Let's see how we can use these tasks in a workgraph.

with WorkGraph("AddMultiplyWithDefinedTasks") as wg:
    wg.inputs = dict.fromkeys(["x", "y", "z"])
    add_outputs = add(x=wg.inputs.x, y=wg.inputs.y)
    mul_outputs = multiply(x=add_outputs.result, y=wg.inputs.z)
    wg.outputs.result = mul_outputs.result

wg.to_html()

# %%
# Using ``add`` and ``multiply`` is quite similar to how functions are used in Python.
# However, unlike the case of ``wg.inputs.x + wg.inputs.y``, calling a ``Task`` returns an outputs namespace.
# This is because in general, a task can have multiple outputs.
# For example, let's compare calling ``add``

add(x=2, y=3)

# %%
# and calling the following:


@task(outputs=["count", "even_numbers"])
def count_even_numbers(numbers):
    even_numbers = [n for n in numbers if n % 2 == 0]
    return {
        "count": len(even_numbers),
        "even_numbers": even_numbers,
    }


count_even_numbers(numbers=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

# %%
# Without explicitly specifying the expected outputs of a task, we get a default ``result`` socket (available output).
# When we do specify our outputs, as done here, they become available at ``<task_name>.outputs.<output_name>`` post-execution.

with WorkGraph("CountEvenNumbers") as wg:
    wg.inputs = dict.fromkeys(["numbers"])
    outputs = count_even_numbers(numbers=wg.inputs.numbers)
    wg.outputs = {
        "count": outputs.count,
        "even_numbers": outputs.even_numbers,
    }

wg.to_html()

wg.submit(
    inputs={
        "numbers": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    },
    wait=True,
)

print("Result:")
print("  Count:", wg.outputs.count.value)
print("  Even numbers:", wg.outputs.even_numbers.value)

# %%
# Great! We can now define functional tasks to cover any Python logic we can think of.
# Let's now think of higher-order workflows, i.e. ones that reuse other workflows as tasks.

# %%
# ``WorkGraph`` as a ``Task``
# ---------------------------
#
# We can expose a ``WorkGraph`` as a ``Task`` by using the ``@task.graph`` decorator.
# This is a great way to encapsulate a workflow and reuse it as a sub-workflow (task) in other workflows.


@task
def random_number_generator(minimum, maximum):
    import random

    return random.randint(minimum, maximum)


@task
def generate_numbers(n):
    return [*range(n)]


@task
def sum_numbers(numbers):
    return sum(numbers)


@task.graph
def SumEvenNumbers(numbers):
    even_numbers_output = count_even_numbers(numbers=numbers)
    return sum_numbers(numbers=even_numbers_output.even_numbers).result


with WorkGraph() as wg:
    wg.inputs = dict.fromkeys(["minimum", "maximum"])
    n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
    numbers = generate_numbers(n=n).result
    wg.outputs = SumEvenNumbers(numbers=numbers).result

wg.to_html()

# %%
# Here we see ``SumEvenNumbers`` as a black-box task.
# To inspect its internal tasks, we can use the ``get_graph`` method as follows:

SumEvenNumbers.get_graph(
    numbers=[*range(42)]
).to_html()  # the input doesn't matter here

# %%
# .. note::
#
#    When using the AiiDA GUI to inspect executed workflows, sub-workflows can be expanded interactively to inspect their internal tasks and connections.

# %%
# .. note::
#
#    Though not strictly required, we name our graph tasks using camel case to distinguish them from regular tasks.
#    If not explicitly overridden, the name of the decorated function will be used as the name of the task when inspecting processes using, for example, ``verdi process list``.

# %%
# Complex workflows
# -----------------
#
# So far, our workflows have been simple and linear.
# However, often we need to create more complex workflows that involve branching, loops, and other control structures.
# Let's explore such a scenario.
# Suppose we want to conditionally sum a set of numbers based on some random input.
# If the input is greather than 50, we filter out the even numbers and sum them.
# Otherwise, we just sum all the numbers.
# Let's try to define this workflow using what we've learned above.

# .. code-block:: python
#
#    with WorkGraph("ComplexWorkflow") as wg:
#        wg.inputs = dict.fromkeys(["minimum", "maximum"])
#        n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
#        numbers = generate_numbers(n=n).result
#        if n < 50:
#            the_sum = sum_numbers(numbers=numbers).result
#        else:
#            the_sum = SumEvenNumbers(numbers=numbers).result
#        wg.outputs.sum = the_sum

# %%
# This seems straightforward and Pythonic. However, if won't work as expected.
# The reason is that the ``if`` statement requires the value of ``n``, which is not available at the time of the workflow definition.
# To handle this conditional logic, and any other flow control for that matter, we must queue it as a graph task.
# Another way to think of it is that any Pythonic control flow statement must be deferred to runtime.
# We do this via task encapsulation.
# Let's do this for the above example.


@task.graph
def ConditionalSum(numbers):
    if len(numbers) < 50:
        return sum_numbers(numbers=numbers).result
    return SumEvenNumbers(numbers=numbers).result


with WorkGraph("ComplexWorkflow") as wg:
    wg.inputs = dict.fromkeys(["minimum", "maximum"])
    n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
    numbers = generate_numbers(n=n).result
    wg.outputs.sum = ConditionalSum(numbers=numbers).result

wg.to_html()

# %%
#
# By treating the conditional logic as a task, there is no need to know the value of ``n`` at workflow definition time, as the conditional statement will only evaluate when the workflow is executed.
# This allows us to define dynamic workflows that decide their flow at runtime.
# To see this more clearly, we can inspect the graph of the ``ConditionalSum`` task given inputs on both ends of the condition.

ConditionalSum.get_graph(numbers=[*range(40)]).to_html()

# %%
ConditionalSum.get_graph(numbers=[*range(60)]).to_html()

# %%
# You can learn more about ``WorkGraph`` flow control in the :doc:`autogen/control-flow` section.
# You can also learn more about the ``@task.graph`` decorator in the :doc:`autogen/graph_builder_concept` section.

# %%
# What’s Next
# -----------
#
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Concepts <../concept/index.rst>`__                     | A brief introduction of ``WorkGraph``’s main concepts. |
# |                                                         |                                                        |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Tutorials <../tutorial/index.rst>`__                   | Real-world examples in computational materials         |
# |                                                         | science and more.                                      |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `HowTo <../howto/index.rst>`__                          | Advanced topics and tips, e.g flow control using       |
# |                                                         | ``if``, ``for``, ``while`` and ``context``.            |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
# | `Migration <migration_from_aiida_core/index.rst>`       | Migration guide for users interested in transferring   |
# |                                                         | aiida-core processes (e.g. ``WorkChain``, ``CalcJob``) |
# |                                                         | to ``WorkGraph``.                                      |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
#
# .. _AiiDA documentation: https://aiida.readthedocs.io/en/stable/
