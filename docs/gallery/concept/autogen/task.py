"""
Task
====

Task is the basic building block of the WorkGraph. A task has inputs,
outputs, and the executor. A task executor can be a ``calcfunction``,
``workfunction``, ``calcjob``, ``Workchain`` or any other Python
function, or even an AiiDA ``ProcessBuilder``.
A task can be created in three ways.

Decorator
---------

Decorate any Python function using the ``task`` decorator. To use the
power of AiiDA (e.g. save the results to a database, keep provenance),
one can use the ``task.calcfunction`` decorator (note that this will,
however, require that the inputs and outputs of your function have to
be instances of ``orm.Node``.

"""

from aiida_workgraph import task
from aiida import orm

# define add task
@task  # this is equivalent to passing no arguments @task()
def add(x, y):
    return x + y


# define multiply calcfunction task
@task.calcfunction  # this is equivalent to passing no arguments @task.calculation()
def multiply(x, y):
    return orm.Float(x + y)


# export the task to html file so that it can be visualized in a browser
add.task().to_html()

# visualize the task in jupyter-notebook
# add.task()


######################################################################
# The input ports (also named sockets) are generated automatically based
# on the function arguments. The default name of the output port is
# ``result``. There are also some built-in ports, like ``_wait`` and
# ``_outputs``. One can create a task instance and inspect its inputs and
# outputs:
#

add1 = add.task()
print("Inputs:", add1.get_input_names())
print("Outputs:", add1.get_output_names())


######################################################################
# If you want to change the name of the output ports, or if there are more
# than one output. You can define the outputs explicitly. For example:


# define the outputs explicitly
@task(outputs=["sum", "diff"])
def add_minus(x, y):
    return {"sum": x + y, "difference": x - y}


print("Inputs:", add_minus.task().get_input_names())
print("Outputs:", add_minus.task().get_output_names())

######################################################################
# One can also add an ``identifier`` to indicates the data type. The data
# type tell the code how to display the port in the GUI, validate the data,
# and serialize data into database.
# We use ``workgraph.Any`` for any data type. For the moment, the data validation is
# experimentally supported, and the GUI display is not implemented. Thus,
# I suggest you to always ``workgraph.Any`` for the port.
#

# define the outputs with identifier
@task(
    outputs=[
        {"name": "sum", "identifier": "workgraph.Any"},
        {"name": "diff", "identifier": "workgraph.Any"},
    ]
)
def add_minus(x, y):
    return {"sum": x + y, "difference": x - y}


######################################################################
# Then, one can use the task inside the WorkGraph:
#

from aiida_workgraph import WorkGraph

wg = WorkGraph()
add_minus1 = wg.add_task(add_minus, name="add_minus1")
multiply1 = wg.add_task(multiply, name="multiply1")
wg.add_link(add_minus1.outputs.sum, multiply1.inputs.x)


######################################################################
# Build from Callable
# -------------------
#
# One can build a task from an already existing Python function.
#

from aiida_workgraph import WorkGraph, build_task

from scipy.linalg import norm

NormTask = build_task(norm)

wg = WorkGraph()
norm_task = wg.add_task(NormTask, name="norm1")
norm_task.to_html()


######################################################################
# The inputs and outputs of the task are automatically generated. One can
# also define the outputs explicitly.
#

NormTask = build_task(norm, outputs=[{"name": "norm", "identifier": "workgraph.Any"}])
wg = WorkGraph()
norm_task = wg.add_task(NormTask, name="norm1")

print("Inputs: ", norm_task.inputs)
print("Outputs: ", norm_task.outputs)

######################################################################
# For specifying the outputs, the most explicit way is to provide a list of dictionaries, as shown above. In addition,
# as a shortcut, it is also possible to pass a list of strings. In that case, WorkGraph will internally convert the list
# of strings into a list of dictionaries in which case, each ``name`` key will be assigned each passed string value.
# Furthermore, also a mixed list of string and dict elements can be passed, which can be useful in cases where multiple
# outputs should be specified, but more detailed properties are only required for some of the outputs. The above also
# applies for the ``outputs`` argument of the ``@task`` decorator introduced earlier, as well as the ``inputs``, given
# that they are explicitly specified rather than derived from the signature of the ``Callable``.  Finally, all lines
# below are valid specifiers for the ``outputs`` of the ``build_task`:
#

NormTask = build_task(norm, outputs=["norm"])
NormTask = build_task(norm, outputs=["norm", "norm2"])
NormTask = build_task(
    norm, outputs=["norm", {"name": "norm2", "identifier": "workgraph.Any"}]
)
NormTask = build_task(
    norm,
    outputs=[
        {"name": "norm", "identifier": "workgraph.Any"},
        {"name": "norm2", "identifier": "workgraph.Any"},
    ],
)

######################################################################
# One can use these AiiDA component directly in the WorkGraph. The inputs
# and outputs of the task is automatically generated based on the input
# and output port of the AiiDA component. In case of ``calcfunction``, the
# default output is ``result``. If there are more than one output task,
# one need to define the outputs explictily.
#

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

wg = WorkGraph()
add1 = wg.add_task(ArithmeticAddCalculation, name="add1")


######################################################################
# Define a Task
# -------------
#
# Create a task class by inheriting from ``Task`` base class.
#

from aiida_workgraph.task import Task


class MyAdd(Task):

    identifier: str = "MyAdd"
    name = "MyAdd"
    node_type = "calcfunction"
    catalog = "Test"

    _executor = {
        "module_path": "aiida_workgraph.executors.test",
        "callable_name": "add",
    }

    def create_sockets(self):
        self.inputs._clear()
        self.outputs._clear()
        inp = self.add_input("workgraph.Any", "x")
        inp = self.add_input("workgraph.Any", "y")
        self.add_output("workgraph.Any", "sum")


######################################################################
# Then, one can use the task by using its identifier.
#

from aiida_workgraph import WorkGraph

wg = WorkGraph()
add1_task = wg.add_task(MyAdd, name="add1")


######################################################################
# One can also register the task in task pool, and then use its
# ``identifer`` directly.
#
# .. code:: python
#
#    wg.add_task("MyAdd", name="add1")
#
