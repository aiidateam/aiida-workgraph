
Welcome to AiiDA WorkGraph's documentation!
===========================================

Provides the third workflow component: ``WorkGraph``, to design flexible node-based workflows using AiiDA.

.. raw:: html

    <div>
        <object type="text/html" data="_static/first_workflow.html" width="100%" height="400px" allowfullscreen="true"></object>
    </div>


In AiiDA, there are two workflow components: `workfunction` and `WorkChain`. Workfunction is easy to implement but it does not support automatic checkpointing, which is important for long-running calculations. Workchain supports automatic checkpointing but it is difficult to implement and also not as flexible as the `workfunction`. AiiDA-WorkGraph provides the third component: `WorkGraph`. It is easy to implement and supports automatic checkpointing. It is also flexible and can be used to design complex workflows. Check this `blog <blog/workgraph_vs_workchain.ipynb>`_ post for more details.


Visit the `Workgraph Collections repository <https://github.com/superstar54/workgraph-collections>`_ to see demonstrations of how to utilize AiiDA Workgraph for different computational codes.



.. toctree::
   :maxdepth: 1
   :caption: Contents:

   quick_start
   installation
   zero_to_hero
   tutorial/index
   howto/index
   built-in/index
   blog/index
   concept/index
   development/index
   faqs



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
