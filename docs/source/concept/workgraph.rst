.. _workgraph:

.. module:: workgraph


===========================================
WorkGraph
===========================================
The :class:`~aiida_workgraph.workgraph.WorkGraph` object is a collection of tasks and links.

Create and launch workgraph
============================
- Create a empty workgraph:

.. code-block:: python

    from aiida_workgraph import WorkGraph, task
    wg = WorkGraph(name="my_first_workgraph")


Create and use `task`.

.. code:: python

    # define add calcfunction task
    @task.calcfunction()
    def add(x, y):
       return x + y

    add1 = wg.add_task(add, name="add1")
    add2 = wg.add_task(add, name="add2")


- Add link between tasks:

.. code-block:: python

    wg.add_link(add1.outputs["result"], add2.inputs["x"])

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
The tasks will be executed when:

- No input task
- All input tasks finish.


Group outputs
=====================================
One can output the results of the tasks as the output of the WorkGraph.

.. code-block:: python

    wg = WorkGraph("test_workgraph_group_outputs")
    wg.add_task(add, "add1", x=2, y=3)
    wg.group_outputs = [{"name": "sum", "from": "add1.result"}]
    wg.submit(wait=True)
    assert wg.process.outputs.sum.value == 5


List of all Methods
===================

.. autoclass:: aiida_workgraph.workgraph.WorkGraph
   :members:
