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
#    pip install aiida-workgraph[widget]
# 
# Start (or restart) the AiiDA daemon if needed:
# 
# .. code:: console
# 
#    verdi daemon start
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
# outputs and an executor. The executor can be a AiiDA ``calcfunction``,
# ``CalcJob``, ``WorkChain`` or any Python function.
# 

from aiida_workgraph import task

# define add task
@task.calcfunction()
def add(x, y):
    return x + y

# define multiply task
@task.calcfunction()
def multiply(x, y):
    return x*y



######################################################################
# Visualize the task
# ^^^^^^^^^^^^^^^^^^
# 
# If you are running in a Jupiter notebook, you can visualize the task,
# including its inputs and output sockets. You can also export the task to
# html file so that it can be visualized in a browser.
# 


# export the task to html file so that it can be visualized in a browser
add.task().to_html()

# comment out the following line to visualize the task in Jupyter Notebook
# add.task()


######################################################################
# The input sockets are generated automatically based on the function
# arguments. The default name of the output socket is ``result``. There
# are also some built-in sockets, like ``_wait`` and ``_outputs``. In case
# of ``calcfunction``, it also has several built-in sockets, such as
# ``metadata``.
# 


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
wg.add_task(multiply, name="multiply1", x = add_task.outputs["result"])

# export the workgraph to html file so that it can be visualized in a browser
wg.to_html()
# comment out the following line to visualize the workgraph in jupyter-notebook
wg


######################################################################
# Submit and view results
# ~~~~~~~~~~~~~~~~~~~~~~~
# 

from aiida.orm import Int

#------------------------- Submit the calculation -------------------
wg.submit(inputs = {"add1": {"x": Int(2),
                             "y": Int(3)
                             },
                    "multiply1": {"y": Int(4)}
                    },
          wait=True)

print("State of WorkGraph:   {}".format(wg.state))
print('Result of add      : {}'.format(wg.tasks["add1"].outputs[0].value))
print('Result of multiply : {}'.format(wg.tasks["multiply1"].outputs[0].value))


######################################################################
# One can also generate the node graph from the AiiDA process:
# 

from aiida_workgraph.utils import generate_node_graph
generate_node_graph(wg.pk)


######################################################################
# Remote job
# ----------
# 
# The ``PythonJob`` is a built-in task that allows users to run Python
# functions on a remote computer.
# 
# In this case, we use define the task using normal function instead of
# ``calcfunction``. Thus, user does not need to install AiiDA on the
# remote computer.
# 

from aiida_workgraph import WorkGraph, task

# define add task
@task()
def add(x, y):
    return x + y

# define multiply task
@task()
def multiply(x, y):
    return x*y

wg = WorkGraph("second_workflow")
# You might need to adapt the code_label to python3 if you use this as your default python
wg.add_task("PythonJob", function=add, name="add", code_label="python")
wg.add_task("PythonJob", function=multiply, name="multiply", x=wg.tasks["add"].outputs[0], code_label="python")

# export the workgraph to html file so that it can be visualized in a browser
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg


######################################################################
# Submit the workgraph
# ~~~~~~~~~~~~~~~~~~~~
# 
# **Code**: We can set the ``computer`` to the remote computer where we
# want to run the job. This will create a code ``python@computer`` if not
# exists. Of course, you can also set the ``code`` directly if you have
# already created the code.
# 
# **Data**: Users can (and is recoomaneded) use normal Python data as
# input. The workgraph will transfer the data to AiiDA data
# (``GeneralData``) using pickle.
# 
# **Python Version**: since pickle is used to store and load data, the
# Python version on the remote computer should match the one used in the
# localhost. One can use conda to create a virtual environment with the
# same Python version. Then activate the environment before running the
# script.
# 
# .. code:: python
# 
#    # For real applications, one can pass metadata to the scheduler to activate the conda environment
#    metadata = {
#        "options": {
#            'custom_scheduler_commands' : 'module load anaconda\nconda activate py3.11\n',
#        }
#    }
# 

from aiida_workgraph.utils import generate_node_graph

#------------------------- Submit the calculation -------------------
# For real applications, one can pass metadata to the scheduler to activate the conda environment
metadata = {
    "options": {
        # 'custom_scheduler_commands' : 'module load anaconda\nconda activate py3.11\n',
        'custom_scheduler_commands' : '',
    }
}

wg.submit(inputs = {"add": {"x": 2, "y": 3,
                            "computer": "localhost",
                            "metadata": metadata},
                    "multiply": {"y": 4,
                            "computer": "localhost",
                            "metadata": metadata}},
          wait=True)
#------------------------- Print the output -------------------------
print("\nResult of multiply is {} \n\n".format(wg.tasks["multiply"].outputs['result'].value))
#------------------------- Generate node graph -------------------
generate_node_graph(wg.pk)


######################################################################
# Use parent folder
# ~~~~~~~~~~~~~~~~~
# 
# The parent_folder parameter allows a task to access the output files of
# a parent task. This feature is particularly useful when you want to
# reuse data generated by a previous computation in subsequent
# computations. In the following example, the multiply task uses the
# ``result.txt`` file created by the add task.
# 
# By default, the content of the parent folder is symlinked to the working
# directory. In the function, you can access the parent folder using the
# relative path. For example, ``./parent_folder/result.txt``.
# 

from aiida_workgraph import WorkGraph, task

# define add task
@task()
def add(x, y):
    z = x + y
    with open("result.txt", "w") as f:
        f.write(str(z))

# define multiply task
@task()
def multiply(x, y):
    with open("parent_folder/result.txt", "r") as f:
        z = int(f.read())
    return x*y + z

wg = WorkGraph("third_workflow")
# You might need to adapt the code_label to python3 if you use this as your default python
wg.add_task("PythonJob", function=add, name="add", code_label="python")
wg.add_task("PythonJob", function=multiply, name="multiply",
             parent_folder=wg.tasks["add"].outputs["remote_folder"],
             code_label="python"
             )

wg.to_html()



######################################################################
# Submit the calculation
# 

#------------------------- Submit the calculation -------------------
wg.submit(inputs = {"add": {"x": 2, "y": 3, "computer": "localhost"},
                    "multiply": {"x": 3, "y": 4, "computer": "localhost"}},
          wait=True)
print("\nResult of multiply is {} \n\n".format(wg.tasks["multiply"].outputs['result'].value))


######################################################################
# CalcJob and WorkChain
# ---------------------
# 
# AiiDA also provide builtin ``CalcJob`` to run a calculation on a remote
# computer. AiiDA community also provides a lot of well-written
# ``calcfunction`` and ``WorkChain``. One can use these AiiDA component
# direclty in the WorkGraph. The inputs and outputs of the task is
# automatically generated based on the input and output port of the AiiDA
# component.
# 
# Here is an example of using the ``ArithmeticAddCalculation`` Calcjob
# inside the workgraph.
#

from aiida_workgraph import WorkGraph
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida.orm import Int, InstalledCode, load_computer



code = InstalledCode(
    computer=load_computer('localhost'),
    filepath_executable='/bin/bash',
    label='add',
    default_calc_job_plugin='core.arithmetic.add',
).store()

wg = WorkGraph("test_add_multiply")
add1 = wg.add_task(ArithmeticAddCalculation, name="add1", x=Int(2), y=Int(3), code=code)
add2 = wg.add_task(ArithmeticAddCalculation, name="add2", y=Int(3), code=code)
wg.add_link(wg.tasks["add1"].outputs["sum"], wg.tasks["add2"].inputs["x"])
wg.to_html()



######################################################################
# Submit the workgraph and wait for the result.
# 

wg.submit(wait=True)
print('Result of task add1: {}'.format(wg.tasks["add2"].outputs["sum"].value))

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


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
@task.calcfunction()
def add(x, y):
    return x + y

# define multiply task
@task.calcfunction()
def multiply(x, y):
    return x*y


# use task.graph_builder decorator, expose the "result" output of "multiply" task
# as the "multiply" output of the `WorkGraph`.
@task.graph_builder(outputs = [{"name": "multiply", "from": "multiply.result"}])
def add_multiply(x, y, z):
    # Create a WorkGraph
    wg = WorkGraph()
    wg.add_task(add, name="add", x=x, y=y)
    wg.add_task(multiply, name="multiply", x=z)
    wg.add_link(wg.tasks["add"].outputs["result"], wg.tasks["multiply"].inputs["y"])
    # don't forget to return the `wg`
    return wg


######################################################################
# Use this graph builder inside a ``WorkGraph``:
# 


from aiida_workgraph import WorkGraph
from aiida.orm import Int

wg = WorkGraph("test_graph_build")
# create a task using the graph builder
add_multiply1 = wg.add_task(add_multiply, x=Int(2), y=Int(3), z=Int(4))
add_multiply2 = wg.add_task(add_multiply, x=Int(2), y=Int(3))
# link the output of a task to the input of another task
wg.add_link(add_multiply1.outputs["multiply"], add_multiply2.inputs["z"])
wg.submit(wait=True)
print("WorkGraph state: ", wg.state)


######################################################################
# Get the result of the tasks:
# 

print('Result of task add_multiply1: {}'.format(add_multiply1.outputs["multiply"].value))

generate_node_graph(wg.pk)


######################################################################
# Start the web server
# ~~~~~~~~~~~~~~~~~~~~
# 
# WorkGraph also provides a web GUI, where you can view and manage the
# workgraph. Open a terminal, and run:
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
# +---------------+------------------------------------------------------+
# | `Conce        | A brief introduction of WorkGraph’s main concepts.   |
# | pts <concept/ |                                                      |
# | index.rst>`__ |                                                      |
# +---------------+------------------------------------------------------+
# | `Tutoria      | Real-world examples in computational materials       |
# | ls <tutorial/ | science and more.                                    |
# | index.rst>`__ |                                                      |
# +---------------+------------------------------------------------------+
# | `             | Advanced topics and tips, e.g flow control using     |
# | HowTo <howto/ | ``if``, ``for``, ``while`` and ``context``.          |
# | index.rst>`__ |                                                      |
# +---------------+------------------------------------------------------+
# 
