.. _workgraph:

.. module:: workgraph


===========================================
WorkGraph
===========================================
The :class:`~aiida_workgraph.workgraph.WorkGraph` object is a collection of nodes and links.

Create and launch workgraph
============================
- Create a empty workgraph:

.. code-block:: python

    from aiida_workgraph import WorkGraph, node
    wg = WorkGraph(name="my_first_workgraph")


Create and use `node`.

.. code:: python

    # define add calcfunction node
    @node.calcfunction()
    def add(x, y):
       return x + y

    add1 = wg.nodes.new(add, name="add1")
    add2 = wg.nodes.new(add, name="add2")


- Add link between nodes:

.. code-block:: python

    wg.links.new(add1.outputs[0], add2.inputs[0])

- Submit the workgraph:

.. code-block:: python

    wg.submit()

Load workgraph from the AiiDA process
=====================================
WorkGraph save its data as a extra attribute into its process, so that one can rebuild the WorkGraph from the process.


.. code-block:: python

    from aiida_workgraph import WorkGraph
    # pk is the process id of a WorkGraph
    WorkGraph.load(pk)

Execute order
===============
The nodes will be executed when:

- No input node
- All input nodes finish.




List of all Methods
===================

.. autoclass:: aiida_workgraph.workgraph.WorkGraph
   :members:
