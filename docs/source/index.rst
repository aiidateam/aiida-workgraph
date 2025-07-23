
AiiDA WorkGraph
===========================================

Efficiently design and manage flexible workflows with AiiDA, featuring an interactive GUI, checkpoints, provenance tracking, and remote execution capabilities.

.. raw:: html

    <div>
        <object type="text/html" data="_static/first_workflow.html" width="100%" height="400px" allowfullscreen="true"></object>
    </div>

Key Features
------------

- **Easy to use**: Create workflows by linking the input and output socket of different tasks.
- **Flexible**: Extend (modify) the workflow by adding (editing) tasks and links, or combine multiple workflows together.
- **Interactive GUI**: Visualize and interact with the workflow using the GUI.
- **Checkpoints**: Save the workflow state, and resume the workflow from the last checkpoint.
- **Provenance**: Track the provenance of the workflow.
- **Remote execution**: Execute the task (Python function, Shell command) on a remote machine.

Check this `blog <blog/workgraph_vs_workchain.ipynb>`_ post for the comparison between WorkGraph and AiiDA's WorkChain.



Sections
========

   .. container:: tocdescr

      .. container:: descr

         :doc:`/autogen/quick_start`
            A quick start guide to get you up and running with AiiDA WorkGraph.

      .. container:: descr

         :doc:`/installation`
            Installation instructions for AiiDA WorkGraph.

      .. container:: descr

         :doc:`/tutorial/index`
            A step-by-step guide to creating a real-world workflow using AiiDA WorkGraph.

      .. container:: descr

         :doc:`/howto/index`
            How-to guides for AiiDA WorkGraph.

      .. container:: descr

         :doc:`/migration_from_aiida_core/index`
            Transition your existing workflows from AiiDA Core to the AiiDA Workgraph

      .. container:: descr

         :doc:`concept/index`
            Concepts and terminologies used in AiiDA WorkGraph.

      .. container:: descr

         :doc:`gui/index`
            Interactive GUI and job menagement of WorkGraph.

      .. container:: descr

         :doc:`development/index`
            Development guide for AiiDA WorkGraph.

      .. container:: descr

         :doc:`gallery`
            Gallery of workflows created using AiiDA WorkGraph.

      .. container:: descr

         :doc:`blog/index`
            Blog posts related to AiiDA WorkGraph.


      .. container:: descr

         :doc:`faqs`
            Frequently asked questions about AiiDA WorkGraph.



.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :hidden:

   autogen/quick_start
   installation
   tutorial/index
   howto/index
   migration_from_aiida_core/index
   concept/index
   gui/index
   development/index
   gallery
   blog/index
   faqs




Indices and tables
--------------------

* :ref:`genindex`
* :ref:`modindex`
