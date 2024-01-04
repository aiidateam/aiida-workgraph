.. _worktree:

.. module:: worktree


===========================================
WorkTree
===========================================
The :class:`~aiida_worktree.worktree.WorkTree` object is a collection of nodes and links.

Create and launch worktree
============================
- Create a empty worktree:

.. code-block:: python

    from aiida_worktree import WorkTree, node
    wt = WorkTree(name="my_first_worktree")


Create and use `node`.

.. code:: python

    # define add calcfunction node
    @node.calcfunction()
    def add(x, y):
       return x + y

    add1 = wt.nodes.new(add, name="add1")
    add2 = wt.nodes.new(add, name="add2")


- Add link between nodes:

.. code-block:: python

    wt.links.new(add1.outputs[0], add2.inputs[0])

- Submit the worktree:

.. code-block:: python

    wt.submit()

Load worktree from the AiiDA process
=====================================
WorkTree save its data as a extra attribute into its process, so that one can rebuild the WorkTree from the process.


.. code-block:: python

    from aiida_worktree import WorkTree
    # pk is the process id of a WorkTree
    WorkTree.load(pk)

Execute order
===============
The nodes will be executed when:

- No input node
- All input nodes finish.




List of all Methods
===================

.. autoclass:: aiida_worktree.worktree.WorkTree
   :members:
