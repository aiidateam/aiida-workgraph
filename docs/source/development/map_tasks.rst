====================
Map task
====================

Copy subgraph
--------------------

Should not deepcopy:

- input value (especially the AiiDA node)


Should update:
- input value from the map
- children
- parent
- input_links


Set context
--------------------
For mapped tasks, should use dict to hold the results


Collect results
--------------------

Inside while zone
--------------------

If the while zone reset, how to handle the copied subgraph, delete or re-use? Note, the node in the subgraph may have been modified by the while zone, so the subgraph may not be the same as the original one.

Edge Cases
--------------------

- Empty subgraph under a MAP node (no children). Just mark the node FINISHED.
- Empty source (no items). Possibly do no clones and mark FINISHED.


Partial failures
--------------------
If one iteration fails, skip other iterations or keep going?
