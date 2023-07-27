.. _node_concept:

===========================================
Node
===========================================

One can create a node in three ways.

Decorator
---------
Decorate a `calcfunction` or `workfunction` with the `@node` decorator.

.. code:: python

   from aiida_worktree import node
   from aiida.engine import calcfunction

   # define add calcfunction node
   @node(outputs = [["General", "sum"]])
   @calcfunction
   def add(x, y):
      return x + y

Then, one can use the node by using its identifier.

.. code:: python

   nt.nodes.new(add.identifier)

Register
--------
Register a already existing `calcfunction`,  `workfunction`, `calcjob`, `Workchain` with the `register_node` function.

.. code:: python

   from aiida_worktree import register_node
   ndata = {"path": "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"}
   add = register_node(ndata)
   # use the node
   nt.nodes.new(add.identifier)


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

   nt.nodes.new("MyAdd")
