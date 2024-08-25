"""
Task
====

Task is the basic building block of the WorkGraph. A task has inputs,
outputs, and the executor. A task executor can be a ``calcfunction``,
``workfunction``, ``calcjob``, ``Workchain`` or any other Python
function. A task can be created in three ways.

Decorator
---------

Decorate any Python function using the ``task`` decorator. To use the
power of AiiDA (e.g. save the results to a database, keep provenance),
one can use the ``task.calcfunction`` decorator.

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
print("Inputs:", add1.inputs.keys())
print("Outputs:", add1.outputs.keys())


######################################################################
# If you want to change the name of the output ports, or if there are more
# than one output. You can define the outputs explicitly. For example:
# ``{"name": "sum", "identifier": "workgraph.Any"}``, where the ``identifier``
# indicates the data type. The data type tell the code how to display the
# port in the GUI, validate the data, and serialize data into database. We
# use ``workgraph.Any`` for any data type. For the moment, the data validation is
# experimentally supported, and the GUI display is not implemented. Thus,
# I suggest you to always ``workgraph.Any`` for the port.
#

# define add calcfunction task
@task(
    outputs=[
        {"name": "sum", "identifier": "workgraph.Any"},
        {"name": "diff", "identifier": "workgraph.Any"},
    ]
)
def add_minus(x, y):
    return {"sum": x + y, "difference": x - y}


print("Inputs:", add_minus.task().inputs.keys())
print("Outputs:", add_minus.task().outputs.keys())


######################################################################
# Then, one can use the task inside the WorkGraph:
#

from aiida_workgraph import WorkGraph

wg = WorkGraph()
add_minus1 = wg.add_task(add_minus, name="add_minus1")
multiply1 = wg.add_task(multiply, name="multiply1")
wg.add_link(add_minus1.outputs["sum"], multiply1.inputs["x"])


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

print("Inputs:")
for input in norm_task.inputs:
    if "." not in input.name:
        print(f"  - {input.name}")
print("Outputs:")
for output in norm_task.outputs:
    if "." not in output.name:
        print(f"  - {output.name}")


######################################################################
# One can use these AiiDA component direclty in the WorkGraph. The inputs
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
    kwargs = ["x", "y"]

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("workgraph.Any", "x")
        inp.add_property("workgraph.Any", "x", default=0.0)
        inp = self.inputs.new("workgraph.Any", "y")
        inp.add_property("workgraph.Any", "y", default=0.0)
        self.outputs.new("workgraph.Any", "sum")

    def get_executor(self):
        return {
            "path": "aiida_workgraph.executors.test",
            "name": "add",
        }


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
