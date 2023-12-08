.. _node_concept:

===========================================
Node
===========================================
Node is the basic building block of the WorkTree. A node has inputs, outputs, and the executor. A node executor can be a `calcfunction`, `workfunction`, `calcjob`, `Workchain` or any other Python function. A node can be created in three ways.

Decorator
---------
Decorate any Python function using the `node` decorator. To use the power of AiiDA (e.g. save result to a database, keep provenance), one can use the AiiDA `calcfunction` decorator with the `node` decorator.

.. code:: python

   from aiida_worktree import node
   from aiida.engine import calcfunction

   # define add calcfunction node
   @node.calcfunction()
   def add(x, y):
      return x + y

Then, one can use the node by using.

.. code:: python

   wt = WorkTree()
   wt.nodes.new(add)

Register
--------
Register a already existing AiiDA `calcfunction`,  `workfunction`, `calcjob`, `Workchain` with the `build_node` function.

.. code:: python

   from aiida_worktree import build_node, WorkTree
   ndata = {"path": "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"}
   add = build_node(ndata)
   # use the node
   wt = WorkTree()
   wt.nodes.new(add)

The inputs and outputs of the node is automatically generated based on the AiiDA process. One can create a node instance and inpsect its inputs and outputs:

.. code:: python

   node = add()
   print("Inputs:")
   for input in node.inputs:
      if "." not in input.name:
         print(f"  - {input.name}")
   print("Outputs:")
   for output in node.outputs:
      if "." not in output.name:
         print(f"  - {output.name}")

Define a Node
-------------
Create a node class by inheriting from `Node` base class.

.. code:: python

   class MyAdd(Node):

      identifier: str = "MyAdd"
      name = "MyAdd"
      node_type = "calcfunction"
      catalog = "Test"
      kwargs = ["x", "y"]

      def create_sockets(self):
          self.inputs.clear()
          self.outputs.clear()
          inp = self.inputs.new("Float", "x")
          inp.add_property("AiiDAFloat", "x", default=0.0)
          inp = self.inputs.new("Float", "y")
          inp.add_property("AiiDAFloat", "y", default=0.0)
          self.outputs.new("Float", "sum")

      def get_executor(self):
          return {
              "path": "aiida_worktree.test",
              "name": "add",
          }

Then, one can use the node by using its identifier.

.. code:: python

   wt.nodes.new("MyAdd")
