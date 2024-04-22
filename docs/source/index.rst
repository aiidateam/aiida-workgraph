
Welcome to AiiDA WorkGraph's documentation!
===========================================

Provides the third workflow component: ``WorkGraph``, to design flexible node-based workflows using AiiDA.

In AiiDA, there are two workflow components: `workfunction` and `WorkChain`. Workfunction is easy to implement but it does not support automatic checkpointing, which is important for long-running calculations. Workchain supports automatic checkpointing but it is difficult to implement and also not as flexible as the `workfunction`. AiiDA-WorkGraph provides the third component: `WorkGraph`. It is easy to implement and supports automatic checkpointing. It is also flexible and can be used to design complex workflows.


Here is a detailed comparison between the ``WorkGraph`` with two AiiDA built-in workflow components. Check this `blog <blog/workgraph_vs_workchain.ipynb>`_ post for more details.


+--------------------------+------------------------+-------------------------------+------------------------+
| Aspect                   | WorkFunction           | WorkChain                     | WorkGraph              |
+==========================+========================+===============================+========================+
| Use Case                 | Short-running          | Long-running                  | Long-running           |
|                          | jobs                   | jobs                          | jobs                   |
+--------------------------+------------------------+-------------------------------+------------------------+
| Checkpointing            | ``No``                 | Yes                           | Yes                    |
+--------------------------+------------------------+-------------------------------+------------------------+
| Execution order          | ``Sequential``         | ``Hybrid Sequential-Parallel``| Directed Acyclic Graph |
+--------------------------+------------------------+-------------------------------+------------------------+
| Non-blocking             | ``No``                 | Yes                           | Yes                    |
+--------------------------+------------------------+-------------------------------+------------------------+
| Implementation           | Easy                   | ``Difficult``                 | Easy                   |
+--------------------------+------------------------+-------------------------------+------------------------+
| Dynamic                  | ``No``                 | ``No``                        | Yes                    |
+--------------------------+------------------------+-------------------------------+------------------------+
| Ready to Use             | Yes                    | ``Need PYTHONPATH``           | Yes                    |
+--------------------------+------------------------+-------------------------------+------------------------+
| Subprocesses Handling    | ``No``                 | Launches & waits              | Launches & waits       |
+--------------------------+------------------------+-------------------------------+------------------------+
| Flow Control             | All                    | `if`, `while`                 | `if`, `while`, `match` |
+--------------------------+------------------------+-------------------------------+------------------------+
| Termination              | ``Hard exit``          | ExitCode                      | ExitCode               |
+--------------------------+------------------------+-------------------------------+------------------------+
| Data Passing             | Direct passing         | Context                       | Link & Context         |
+--------------------------+------------------------+-------------------------------+------------------------+
| Output Recording         | Limited support        | Out & validates               | Out                    |
+--------------------------+------------------------+-------------------------------+------------------------+
| Port Exposing            | Limited support        | Manual & automatic            | Manual                 |
+--------------------------+------------------------+-------------------------------+------------------------+

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   quick_start
   installation
   tutorial/index
   howto/index
   blog/index
   concept/index
   faqs



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
