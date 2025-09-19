"""
Task
====

Task is the basic building block of the WorkGraph. A task has inputs,
outputs, and the executor. A task executor can be a Python function, AiiDA components
(``calcfunction``, ``workfunction``, ``calcjob``, ``Workchain``, and ``ProcessBuilder``).
A task can be created in three ways.

Decorator
---------

Decorate any Python function using the ``task`` decorator.

"""

from aiida_workgraph import task
from aiida import orm


# define add task
@task  # this is equivalent to passing no arguments @task()
def add(x, y):
    return x + y


# define AiiDA calcfunction task
@task.calcfunction  # this is equivalent to passing no arguments @task.calculation()
def multiply(x, y):
    return orm.Float(x + y)


# export the task to html file so that it can be visualized in a browser
add()._node.to_html()

# visualize the task in jupyter-notebook
# add()._node


######################################################################
# The input ports (also named sockets) are generated automatically based
# on the function arguments. The default name of the output port is
# ``result``. There are also some built-in ports, like ``_wait`` and
# ``_outputs``. One can create a task instance and inspect its inputs and
# outputs:
#

add1 = add()._node
print('Inputs:', add1.get_input_names())
print('Outputs:', add1.get_output_names())


######################################################################
# If you want to change the name of the output ports, or if there are more
# than one output. You can define the outputs explicitly.
# Please refer to the `Socket <./socket_concept.rst>`__ for more details.


######################################################################
# Then, one can use the task inside the WorkGraph:
#

from aiida_workgraph import WorkGraph

wg = WorkGraph()
add1 = wg.add_task(add, name='add1')
multiply1 = wg.add_task(multiply, name='multiply1')
wg.add_link(add1.outputs.result, multiply1.inputs.x)


######################################################################
# Build from Callable
# -------------------
#
# One can build a task from an already existing Python function.
#

from aiida_workgraph import WorkGraph, task

from scipy.linalg import norm

NormTask = task()(norm)

wg = WorkGraph()
norm_task = wg.add_task(NormTask, name='norm1')
norm_task.to_html()


######################################################################
# The inputs of the task are automatically generated. However, one need to define the outputs explicitly if there are more than one output.
#


# Define a function with multiple outputs
def calculate_stats(data):
    """
    Calculates the mean and standard deviation of an array.
    """
    import numpy as np

    mean_val = np.mean(data)
    std_val = np.std(data)
    return mean_val, std_val


NormTask = task(outputs=['mean', 'std'])(calculate_stats)
wg = WorkGraph()
norm_task = wg.add_task(NormTask, name='calculate_stats')

print('Inputs: ', norm_task.inputs)
print('Outputs: ', norm_task.outputs)


######################################################################
# One can use AiiDA components directly in the WorkGraph. The inputs
# and outputs of the task is automatically generated based on the input
# and output port of the AiiDA component. In case of ``calcfunction``, the
# default output is ``result``. If there are more than one output task,
# one need to define the outputs explictily.
#

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

wg = WorkGraph()
add1 = wg.add_task(ArithmeticAddCalculation, name='add1')


######################################################################
# Define a Task
# -------------
#
# Create a task class by inheriting from ``Task`` base class.
#

from aiida_workgraph.task import Task


class MyAdd(Task):
    identifier: str = 'MyAdd'
    name = 'MyAdd'
    node_type = 'calcfunction'
    catalog = 'Test'

    _executor = {
        'module_path': 'aiida_workgraph.executors.test',
        'callable_name': 'add',
    }

    def update_sockets(self):
        self.inputs._clear()
        self.outputs._clear()
        _ = self.add_input('workgraph.Any', 'x')
        _ = self.add_input('workgraph.Any', 'y')
        self.add_output('workgraph.Any', 'sum')


######################################################################
# Then, one can use the task by using its identifier.
#

from aiida_workgraph import WorkGraph

wg = WorkGraph()
add1_task = wg.add_task(MyAdd, name='add1')


######################################################################
# One can also register the task in task pool, and then use its
# ``identifer`` directly.
#
# . code:: python
#
#    wg.add_task("MyAdd", name="add1")
#
