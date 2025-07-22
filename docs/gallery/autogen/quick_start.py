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
#
#   - tasks to encapsulate a process
#   - sub-workflows as tasks to encapsulate several sub-processes
#
# - an advanced example demonstrating control flow in ``WorkGraph``
#
# .. note::
#
#    This quick start tutorial assumes basic AiiDA knowledge.
#    To learn more about AiiDA, please visit the `AiiDA documentation`_.

# %%
# Installation
# ------------
#
# Let's first install ``aiida-workgraph`` (and ``aiida-core`` by extension):
#
# .. code:: console
#
#    $ pip install aiida-workgraph

# %%
# Setup
# -----
#
# To interact with the AiiDA database, we need to load an AiiDA profile.
# If you haven't configured one yet, you can do so by running the following command:
#
# .. code:: console
#
#    $ verdi presto
#
# To load your profile, add the following code to your script or Jupyter notebook:

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel("REPORT")
# sphinx_gallery_end_ignore

from aiida import load_profile

load_profile()

# %%
# .. note::
#
#    AiiDA also provides a pre-configured (pre-loaded profile) shell.
#    You can launch it with:
#
#    .. code:: console
#
#       $ verdi shell
#
# .. tip::
#
#    Some features of ``WorkGraph`` are best demonstrated interactively.
#    We recommend using the AiiDA GUI for this.
#    You can learn how to run and use the GUI in the :doc:`../gui/web` section.

# %%
# Simple workflow
# ---------------
#
# Suppose you want to compute ``(x + y) * z`` and record both operations.
# The simplest way to do this in ``WorkGraph`` is as follows:

from aiida_workgraph import WorkGraph

with WorkGraph("AddMultiply") as wg:
    wg.inputs = dict.fromkeys(("x", "y", "z"))
    wg.outputs.result = (wg.inputs.x + wg.inputs.y) * wg.inputs.z

# %%
# We've defined our first workflow using ``WorkGraph`` 🎉
#
# The workflow sets up three input parameters, ``x``, ``y``, and ``z``.
# It then sets up tasks for adding ``x`` and ``y`` (in parentheses), and multiplying the result by ``z``.
# Finally, it assigns the result to the ``result`` output socket.
#
# .. note::
#
#    To learn more about graph-level inputs and output, visit the :doc:`../howto/autogen/graph_level` how-to section.
#
# We can visualize the workgraph to inspect its tasks and their connections.

wg.to_html()

# %%
# Let's run it with some inputs:

wg.run(inputs={"x": 2, "y": 3, "z": 4})

# %%
# .. tip::
#
#    If you've configure a message broker in AiiDA, you can also submit the workgraph
#    (queue it for execution) with
#
#    .. code:: python
#
#       wg.submit(inputs=..., wait=True)
#
#    You can use ``wait=True`` to block the terminal or notebook until the workflow finishes.


# %%
# Let's have a look at the result:

print("Result:", wg.outputs.result.value)

# %%
# We can also access the results of the individual tasks:

print("Result of addition:", wg.tasks.op_add.outputs.result.value)
print("Result of multiplication:", wg.tasks.op_mul.outputs.result.value)

# %%
# .. note::
#
#    Many Python operators are supported in ``WorkGraph`` and are automatically converted to tasks with the name ``op_<operator>``, e.g. ``op_add``, ``op_mul``, ``op_lt``, ``op_eq``, and more.

# %%
# .. tip::
#
#    If you're running the AiiDA GUI, you can visualize the executed workflow interactively.
#    Click on the ``PK`` field of the submitted workflow (look for *WorkGraph<AddMultiply>*) to view its details.
#
# So far so good! We have a simple workflow, can visualize it, run/submit it with inputs, and inspect its outputs.
# However, in practice, tasks may require more advanced logic.
# How do we handle these?

# %%
# The ``@task`` decorator
# -----------------------
#
# We can decorate any Python function with the ``@task`` decorator to turn it into a task.

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


# %%
# We can then call these functional tasks directly in a workgraph as follows:

with WorkGraph("AddMultiplyWithDefinedTasks") as wg:
    wg.inputs = dict.fromkeys(("x", "y", "z"))
    the_sum = add(x=wg.inputs.x, y=wg.inputs.y).result
    the_product = multiply(x=the_sum, y=wg.inputs.z).result
    wg.outputs.result = the_product

wg.to_html()

# %%
# Compare the graph above with the one we defined earlier.
# Other than the task names, they are identical, as to be expected.
#
# Functional tasks are used just as they do in Python.
# The main difference is that instead of values, we pass inputs by reference (sockets).
# This creates a link, or dependency, between input/output sockets.
# Moreover, unlike the case of ``wg.inputs.x + wg.inputs.y``, calling a ``Task`` returns a socket namespace.
# This is because in general, a task can have multiple outputs.
# For example, let's compare calling ``add``:

add(x=2, y=3)

# %%
# with calling the following:


@task(outputs=["count", "even_numbers"])
def count_even_numbers(numbers):
    even_numbers = [n for n in numbers if n % 2 == 0]
    return len(even_numbers), even_numbers


count_even_numbers(numbers=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

# %%
# We can see that without explicitly specifying the expected outputs of a task, we get the default ``result`` socket.
#
# .. note::
#
#    When we do specify our output sockets explicitly, they become available at ``<task_name>.outputs.<output_name>`` post-execution.
#
# .. note::
#
#    To learn more about sockets, see the :doc:`../concept/autogen/socket_concept` concept section.
#
# Great! We can now define functional tasks to cover any Python logic we can think of.
# Let's now consider higher-order workflows.

# %%
# The ``@task.graph`` decorator
# -----------------------------
#
# We can expose a ``WorkGraph`` as a task by using the ``@task.graph`` decorator.
# This is useful for task encapsulation and reusability, as well as for defining dynamic (runtime-input-dependent) workflows.
# Let's see how this works in practice by defining a workflow to sum the even numbers of a number list of random size.


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
    even_numbers = count_even_numbers(numbers=numbers).even_numbers
    return sum_numbers(numbers=even_numbers).result


with WorkGraph("SumRandomEvenNumbers") as wg:
    wg.inputs = dict.fromkeys(("minimum", "maximum"))
    n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
    numbers = generate_numbers(n=n).result
    wg.outputs.sum = SumEvenNumbers(numbers=numbers).result

wg.to_html()

# %%
# Here we see the *SumEvenNumbers* task as a black-box.
# To inspect its internal tasks, we can use the ``build_graph`` method with any input:

SumEvenNumbers.build_graph(numbers=[*range(42)]).to_html()

# %%
# .. tip::
#
#    When using the AiiDA GUI to inspect executed workflows, sub-workflows can be expanded interactively to inspect their internal tasks and connections.
#
# .. note::
#
#    Though not strictly required, we name our graph tasks using camel case to distinguish them from regular tasks.
#    If not explicitly overridden, the name of the decorated function will be used as the name of the task when inspecting processes using, for example, ``verdi process list``.
#
# .. note::
#
#    To learn more about workflow composition, see the :doc:`../howto/autogen/combine_workgraphs` how-to section.
#
# Nice! We can now define reusable workflows, encapsulating any number of tasks in linear fashion.
# But what if our logic isn't linear? What if its iterative, or conditional?

# %%
# Control flow
# ------------
#
# Suppose we want to sum a set of numbers conditionally, depending on how many numbers we have.
# Let's try to define an example workflow using what we've learned so far:
#
# .. code:: python
#
#    with WorkGraph("ComplexWorkflow") as wg:
#        wg.inputs = dict.fromkeys(("minimum", "maximum"))
#        n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
#        numbers = generate_numbers(n=n).result
#        if n < 50:
#            the_sum = sum_numbers(numbers=numbers).result
#        else:
#            the_sum = SumEvenNumbers(numbers=numbers).result
#        wg.outputs.sum = the_sum

# %%
# This seems straightforward.
# However, **it won't work as expected**.
# The reason is that the ``if`` statement requires the value of ``n``, which is not yet available at the time of the workflow definition.
#
# How do we do this then? Simple - we think in tasks!
#
# To handle such conditional logic, or any other flow control for that matter, we must queue it as a dynamic graph task.
# For this, **we need task encapsulation**.
#
# .. important::
#
#    Control flow statements must be deferred to runtime.
#
# Let's do this now for the above example:


@task.graph
def ConditionalSum(numbers):
    if len(numbers) < 50:
        return sum_numbers(numbers=numbers).result
    return SumEvenNumbers(numbers=numbers).result


with WorkGraph("ComplexWorkflow") as wg:
    wg.inputs = dict.fromkeys(("minimum", "maximum"))
    n = random_number_generator(wg.inputs.minimum, wg.inputs.maximum).result
    numbers = generate_numbers(n=n).result
    wg.outputs.sum = ConditionalSum(numbers=numbers).result

wg.to_html()

# %%
#
# By treating the conditional logic as a task, there is no need to know the value of ``n`` at workflow definition time, as the conditional statement will only evaluate when the workflow is executed.
# This allows us to define dynamic workflows that **decide their flow at runtime**.
#
# To see this more clearly, we can inspect the graph of the ``ConditionalSum`` task given inputs on both ends of the condition:

ConditionalSum.build_graph(numbers=[*range(40)]).to_html("html/ConditionalSumTrue.html")

# %%
ConditionalSum.build_graph(numbers=[*range(60)]).to_html(
    "html/ConditionalSumFalse.html"
)

# %%
# .. note::
#
#    To learn more about ``WorkGraph`` flow control, visit the :doc:`../howto/autogen/control-flow` how-to section.
#    To learn more about the ``@task.graph`` decorator, visit the :doc:`../concept/autogen/graph_builder_concept` concept section.


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
# | `Migration <migration_from_aiida_core/index.rst>`__     | Migration guide for users interested in transferring   |
# |                                                         | aiida-core processes (e.g. ``WorkChain``, ``CalcJob``) |
# |                                                         | to ``WorkGraph``.                                      |
# |                                                         |                                                        |
# +---------------------------------------------------------+--------------------------------------------------------+
#
# .. _AiiDA documentation: https://aiida.readthedocs.io/en/stable/

# sphinx_gallery_start_ignore
set_aiida_loglevel("ERROR")
# sphinx_gallery_end_ignore
