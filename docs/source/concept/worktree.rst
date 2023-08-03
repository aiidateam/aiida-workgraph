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

    from aiida_worktree import WorkTree
    wt = WorkTree(name="my_first_worktree")

- Add nodes by using the node identifier.

.. code-block:: python

    float1 = wt.nodes.new("AiiDAFloat", name = "float1")
    float2 = wt.nodes.new("AiiDAFloat", name = "float2")
    sumdiff1 = wt.nodes.new("AiiDASumDiff", name = "sumdiff1")

- Add link between nodes:

.. code-block:: python

    wt.links.new(float1.outputs[0], sumdiff1.inputs[0])
    wt.links.new(float2.outputs[0], sumdiff1.inputs[1])

- Submit the worktree:

.. code-block:: python

    wt.submit()


Execute order
===============
The nodes will be executed when:

- No input node
- All input nodes finish.


List of all Methods
===================

.. autoclass:: aiida_worktree.worktree.WorkTree
   :members:
