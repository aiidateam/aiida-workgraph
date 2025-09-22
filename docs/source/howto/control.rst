.. _topics:control:

*********************************************
How to control (pause, play, kill) the task
*********************************************


Command line
------------
You can control the task from the command line using the `workgraph` command. For example, to pause the task:

.. code-block:: bash

    # pause the task
    workgraph task pause <workgraph_pk> <task_name>
    # play the task
    workgraph task play <workgraph_pk> <task_name>
    # kill the task
    workgraph task kill <workgraph_pk> <task_name>

If you don't know the workgraph id, you can list the running processes:

.. code-block:: bash

    verdi process list

If you don't know the task name, you can list the tasks in the workgraph:

.. code-block:: bash

    workgraph task list <workgraph_pk>

Graphical User Interface
------------------------
.. warning::
   **This feature is experimental.** The API is subject to change in future releases. We welcome your feedback on its functionality.

You can control the task from the graphical user interface. Click on the task and then click on the control buttons. You can pause, play, or kill the task.


.. image:: ../_static/images/control_buttons.png
    :width: 300px
    :align: center
