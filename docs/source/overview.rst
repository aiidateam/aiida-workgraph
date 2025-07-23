========
Overview
========

AiiDA-WorkGraph is a powerful Python library built on the AiiDA framework to streamline the creation, management, and execution of scientific workflows.
It combines a clean Pythonic interface with robust data provenance, high-throughput capabilities, and remote execution—supporting scalable and reproducible research.

Key Features
============

Pythonic workflows
------------------
Define workflows using standard Python functions enhanced with decorators. It's simple and intuitive for new users.

.. code-block:: python

    from aiida_workgraph import task

    @task
    def add(x, y):
        return x + y

    @task
    def multiply(x, y):
        return x * y

    @task.graph()
    def add_multiply(x, y, z):
        sum = add(x, y).result
        product = multiply(sum, z).result
        return product

    add_multiply.run(x=2, y=3, z=4)



Remote execution
----------------
Run Python functions or shell commands on remote machines.

See: `Run calculations remotely <../howto/autogen/remote_job>`_


Shell command example:

.. code-block:: python

    from aiida_workgraph import shelljob

    remote_computer = orm.load_computer("remote_computer_label")

    outputs = shelljob(
        command="date", # the command to run
        metadata={"computer": remote_computer},
    )



Parallel tasks
--------------
Launch multiple tasks in parallel without writing concurrent code.

See: `Run tasks in parallel <../howto/autogen/parallel>`_



.. code-block:: python

    @task.graph()
    def parallel_add(x, N):
        results = {}
        # Launch N parallel tasks to add x with each index
        for i in range(N):
            results[f"add_{i}"] = add(x, i).result
        return results



Provenance tracking
-------------------
Automatically records data and process provenance for full reproducibility and traceability.

Checkpointing
-------------
Save and resume workflow execution from the last checkpoint.

Dynamic execution
-----------------
Adapt workflows at runtime using ``if``, ``while``, and other control structures.

See: `Control flow in WorkGraph <../howto/autogen/control-flow>`_

`if` condition example:

.. code-block:: python

    @task.graph()
    def conditional_workflow(x):
        if x > 0:
            return add(x, 10).result
        else:
            return multiply(x, 10).result


`while` loop example using recursion:

.. code-block:: python

    @task.graph()
    def recursive_workflow(x):
        if x <= 0:
            return x
        else:
            return recursive_workflow(x - 1).result


High-throughput
---------------
Efficiently manage thousands of tasks for large-scale computations.

Error handling
--------------
Recover from failures with built-in retries and error-catching mechanisms.

See: `Write error-resistant workflows <../howto/autogen/error_resistant>`_

Reusable workflows
------------------
Encapsulate and reuse tasks and sub-workflows in larger pipelines.

See: `Combine workgraphs <../howto/autogen/combine_workgraphs>`_

.. code-block:: python

    @task.graph()
    def main_workflow(x, y):
        sum1 = add(x, y).result
        # call the reusable add_multiply workflow
        result1 = add_multiply(sum1, 2, 3).result
        return result1


Interactive GUI
---------------
Visualize and monitor workflows via a user-friendly web interface.

See: `WorkGraph GUI <../gui/autogen/web>`_

.. image:: ./_static/images/web-detail.png

Event-driven logic
------------------
Trigger task execution based on external events for adaptive workflows.
Some possible use cases include:

- **Time-based events**: Start a task at a specified time
- **File-based events**: Execute a task when a particular file exists

Here is an example of defining a monitor task that checks if a certain time has passed:

.. code-block:: python

    @task.monitor
    def time_monitor(time):
        """Monitor a time condition."""
        import datetime
        return datetime.datetime.now() > datetime.datetime.fromisoformat(time.value)


Zone-based control
------------------
Use ``If``, ``While``, and ``For`` zones to explicitly define logic blocks in the workflow graph.




Node-graph editing
------------------
Design workflows by connecting task inputs and outputs like a flowchart.

See: `Node-graph programming <../howto/autogen/node_graph_programming>`_

What's Next?
============

Explore the following resources to get started or dive deeper into AiiDA-WorkGraph:

+-----------------------------------------+------------------------------------------------------+
| `Quick Start <./quick_start.rst>`__     | Get up and running with a simple workflow example.   |
+-----------------------------------------+------------------------------------------------------+
| `Concepts <../concept/index.rst>`__     | Learn the core concepts behind AiiDA-WorkGraph.      |
+-----------------------------------------+------------------------------------------------------+
| `Tutorials <../tutorial/index.rst>`__   | Discover real-world examples in computational        |
|                                         | materials science and other domains.                 |
+-----------------------------------------+------------------------------------------------------+
| `How-To Guides <../howto/index.rst>`__  | Master advanced topics like control flow with        |
|                                         | ``if``, ``for``, ``while``, and ``context``.         |
+-----------------------------------------+------------------------------------------------------+
