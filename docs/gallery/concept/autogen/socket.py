"""
Port (Socket)
=============

In WorkGraph, we use ``sockets`` to indicate the type of data that can
be transferred from one task to another. This is similar to AiiDA’s
``port``. We will use the name ``port`` to reuse the concepts already in
AiiDA as much as possible. Their differences will be introduced later.

Usually, the ports are created automatically from an AiiDA component
(e.g., WorkChain), or generated automatically based on the function
arguments. There are also some built-in ports(sockets), like ``_wait``
and ``_outputs``.

"""

from aiida_workgraph import task
from aiida.manage import load_profile

load_profile()


@task.calcfunction()
def multiply(x, y):
    return x * y


print("Input ports: ", multiply.task().get_input_names())
print("Output ports: ", multiply.task().get_output_names())

multiply.task().to_html()


######################################################################
# If you want to change the name of the output ports, or if there are more
# than one output. You can define the outputs explicitly.
#

from aiida_workgraph import task


@task(
    outputs=[
        {"identifier": "workgraph.Any", "name": "sum"},
        {"identifier": "workgraph.Any", "name": "diff"},
    ]
)
def add_minus(x, y):
    return {"sum": x + y, "difference": x - y}


print("Input ports: ", add_minus.task().get_input_names())
print("Ouput ports: ", add_minus.task().get_output_names())
add_minus.task().to_html()


######################################################################
# Two values are needed to define a port, e.g.,
# ``{"identifier": "General", "name": "sum"}``, where the ``identifier``
# indicates the data type, and the name of the port. We use ``General``
# for any data type.
#
# Assign socket type based on typing hints
# ----------------------------------------
#
# The type hints in the function signature can be used to assign the
# socket type. The following table shows the mapping between the type
# hints and the socket type.
#
# ============= ===============
# Type hint     Socket type
# ============= ===============
# ``orm.Int``   ``AiiDAInt``
# ``orm.Str``   ``AiiDAString``
# ``orm.Float`` ``AiiDAFloat``
# ``orm.Bool``  ``AiiDABool``
# ============= ===============
#

from aiida_workgraph import task


@task.calcfunction()
def add(x: int, y: float) -> float:
    return x + y


print("inputs: ", add.task().inputs)


######################################################################
# Data validation (**Experimental**)
# ----------------------------------
#
# One can use the class of the data directly when defining the port.
#
# **For the moment, data validation is experimentally supported.** Thus, I
# suggest you always use ``workgraph.Any`` for the port.
#

from aiida_workgraph import task
from aiida import orm


@task.calcfunction(
    inputs=[
        {"identifier": orm.Int, "name": "x"},
        {"identifier": orm.Float, "name": "y"},
    ],
    outputs=[{"idenfier": orm.Float, "name": "result"}],
)
def add(x, y):
    result = x + y
    return result


######################################################################
# Advanced concept of Socket
# --------------------------
#
# In the GUI of node graph programming, a socket is displayed as a circle
# only. In order to set the value for a socket directly in the GUI, one
# can add a property to it. A property is the data that can be
# displayed/edited in the GUI directly, which is usually a simple data
# type, such as int, string, boolean, etc.
#
# Property
# ~~~~~~~~
#
# A socket can has a property. The data of the property will be used when
# there is no connection to the input port. The property can be added when
# define a custom port. Or it can be added later by using ``add_property``
# method.
#
# In ``aiida-workgraph``, all socket must have a property. The value of
# the property will be used when there is no connection to the input port.
#


def create_sockets(self):
    # create a General port.
    inp = self.add_input("workgraph.Any", "symbols")
    # add a string property to the port with default value "H".
    inp.add_property("String", "default", default="H")


######################################################################
# Serialization
# ~~~~~~~~~~~~~
#
# If you use non-AiiDA data as inputs/outputs of a ``Normal`` task, the
# data type of the socket will also indicate how to serialize data and
# deserialize the data.
#
