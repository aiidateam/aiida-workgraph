"""
Quick Start
===========

"""


######################################################################
# Introduction
# ------------
#
# To run this tutorial, you need to install ``aiida-workgraph``. Open a
# terminal and run:
#
# .. code:: console
#
#    pip install aiida-workgraph
#
# If you haven't configured an AiiDA profile, you can run the following
#
# .. code:: console
#
#    verdi presto
#
# Load the AiiDA profile.
#

from aiida import load_profile

load_profile()


######################################################################
# First workflow
# --------------
#
# Suppose we want to calculate ``(x + y) * z`` in two steps:
#
# -  add ``x`` and ``y``
# -  then multiply the result with ``z``.
#


######################################################################
# Create task
# ~~~~~~~~~~~
#
# Task is the basic building block of the WorkGraph. A task has inputs,
# outputs and an executor.
#

from aiida_workgraph import task

# define add task
@task()
def add(x, y):
    return x + y


# define multiply task
@task()
def multiply(x, y):
    return x * y


######################################################################
# Create the workflow
# ~~~~~~~~~~~~~~~~~~~
#
# Three steps:
#
# -  create an empty ``WorkGraph``
# -  add tasks: ``add`` and ``multiply``.
# -  link the output of the ``add`` task to the ``x`` input of the
#    ``multiply`` task.
#

from aiida_workgraph import WorkGraph

wg = WorkGraph("add_multiply_workflow")
add_task = wg.add_task(add, name="add1")
# link the output of the `add` task to one of the `x` input of the `multiply` task.
wg.add_task(multiply, name="multiply1", x=add_task.outputs.result)

# export the workgraph to html file so that it can be visualized in a browser
wg.to_html()
# comment out the following line to visualize the workgraph in jupyter-notebook
# wg


######################################################################
# Run and view results
# ~~~~~~~~~~~~~~~~~~~~~~~
#

from aiida.orm import Int

# ------------------------- Run the calculation -------------------
wg.run(
    inputs={"add1": {"x": Int(2), "y": Int(3)}, "multiply1": {"y": Int(4)}},
)

print("State of WorkGraph:   {}".format(wg.state))
print("Result of add      : {}".format(wg.tasks.add1.outputs.result.value))
print("Result of multiply : {}".format(wg.tasks.multiply1.outputs.result.value))


######################################################################
# CalcJob and WorkChain
# ---------------------
#
# One can use AiiDA components (``CalcJob`` and ``WorkChain``)
# direclty in the WorkGraph. The inputs and outputs of the task is
# automatically generated based on the input and output port of the AiiDA
# components.
#
# Here is an example of using the ``ArithmeticAddCalculation`` Calcjob
# inside the workgraph.
#

from aiida_workgraph import WorkGraph
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida.orm import Int, InstalledCode, load_computer, load_code
from aiida.common.exceptions import NotExistent


try:
    code = load_code("add@localhost")  # The computer label can also be omitted here
except NotExistent:
    code = InstalledCode(
        computer=load_computer("localhost"),
        filepath_executable="/bin/bash",
        label="add",
        default_calc_job_plugin="core.arithmetic.add",
    ).store()

wg = WorkGraph("test_add_multiply")
add1 = wg.add_task(ArithmeticAddCalculation, name="add1", x=Int(2), y=Int(3), code=code)
add2 = wg.add_task(ArithmeticAddCalculation, name="add2", y=Int(3), code=code)
wg.add_link(wg.tasks.add1.outputs.sum, wg.tasks.add2.inputs.x)
wg.to_html()


######################################################################
# Run the workgraph and wait for the result.
#

wg.run()
print("\nResult of task add1: {}".format(wg.tasks.add2.outputs.sum.value))


######################################################################
# One can also create task from an AiiDA process builder directly.
#

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

builder = ArithmeticAddCalculation.get_builder()
builder.code = code
builder.x = Int(2)
builder.y = Int(3)

wg = WorkGraph("test_set_inputs_from_builder")
add1 = wg.add_task(builder, name="add1")


######################################################################
# Graph builder
# -------------
#
# A ``WorkGraph`` is a group of tasks. One can treat a ``WorkGraph`` as a
# single task, and expose the inputs and outputs of the ``WorkGraph``.
# This allow you to write:
#
# -  nested workflows
# -  dynamic workflow based on the input values. For example, if you want
#    to use ``if`` and ``for`` to create the tasks, or repeat a
#    calculation until it converges.
#
# The ``Graph Builder`` allow user to create a dynamic workflow based on
# the input value, as well as nested workflows. Here is an example of
# nested workflow:
#

from aiida_workgraph import WorkGraph, task

# define add task
@task()
def add(x, y):
    return x + y


# define multiply task
@task()
def multiply(x, y):
    return x * y


# use task.graph decorator, expose the "result" output of "multiply" task
# as the "multiply" output of the `WorkGraph`.
@task.graph()
def add_multiply(x, y, z):
    outputs1 = add(x=x, y=y)
    outputs2 = multiply(x=z, y=outputs1.result)
    return outputs2.result


######################################################################
# Use this graph builder inside a ``WorkGraph``:
#


from aiida_workgraph import WorkGraph
from aiida.orm import Int

with WorkGraph("test_graph_build") as wg:
    # create a task using the graph builder
    outputs1 = add_multiply(x=Int(2), y=Int(3), z=Int(4))
    # link the output of a task to the input of another task
    outputs2 = add_multiply(x=Int(2), y=Int(3), z=outputs1.result)
    wg.run()


######################################################################
# Get the result of the tasks:
#

print("WorkGraph state: ", wg.state)
print("Result of task add_multiply1: {}".format(outputs2.result.value))


######################################################################
# Start the web server
# ~~~~~~~~~~~~~~~~~~~~
#
# WorkGraph also provides a web GUI, where you can view and manage the
# workgraph. To use the web ui, first install the web ui package:
#
# ::
#
#    pip install aiida-workgraph-web-ui
#
# Open a terminal, and run:
#
# ::
#
#    workgraph web start
#
# Then visit the page http://127.0.0.1:8000/workgraph, you can view the
# workgraph later from here. You should find all the submited workgraph,
# e.g., the ``first_workflow`` WorkGraph. Please click the pk and view the
# workgraph.
#


######################################################################
# What’s Next
# -----------
#
# +-----------------------------------------+------------------------------------------------------+
# | `Concepts <../concept/index.rst>`__     | A brief introduction of WorkGraph’s main concepts.   |
# |                                         |                                                      |
# |                                         |                                                      |
# +-----------------------------------------+------------------------------------------------------+
# | `Tutorials <../tutorial/index.rst>`__   | Real-world examples in computational materials       |
# |                                         | science and more.                                    |
# |                                         |                                                      |
# +-----------------------------------------+------------------------------------------------------+
# | `HowTo <../howto/index.rst>`__          | Advanced topics and tips, e.g flow control using     |
# |                                         | ``if``, ``for``, ``while`` and ``context``.          |
# |                                         |                                                      |
# +-----------------------------------------+------------------------------------------------------+
#
#
