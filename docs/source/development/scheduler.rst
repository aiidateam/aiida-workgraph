====================
Scheduler
====================

Scheduler Design Overview
=========================

Architecture
------------

- The ``workgraph`` scheduler operates with its own dedicated queue, separate from AiiDA’s default process queue.
- Users can submit processes (e.g., ``WorkGraph``, ``CalcJob``) to this scheduler. When submitted:

  1. The process is stored in the AiiDA database as usual.
  2. **However**, no message is sent to AiiDA’s default process queue — meaning the process is not immediately picked up by the daemon.
  3. Instead, a message is sent to the scheduler’s **custom queue**.

Execution Flow
--------------

- The scheduler listens to its queue and decides whether a process should be launched based on the current load and configuration.
- If the conditions are met:

  - It sends a message to the standard AiiDA process queue.
  - The process is then executed by the regular AiiDA daemon runner.

Monitoring and Feedback
-----------------------

- When launching a process, the scheduler registers a **broadcast subscriber**.
- This subscriber is notified when a process finishes, allowing the scheduler to check and potentially launch more waiting processes.

Priority Management
-------------------

- The scheduler determines which process to launch next based on **priority**.
- Each submitted process is assigned a priority, which decreases incrementally with each new submission (i.e., first = 0, next = -1, etc.).
- For a ``WorkGraph`` submitted to the scheduler:

  - All its child processes are also assigned to the same scheduler.
  - These child processes inherit the same priority.

Dynamic Control
---------------

- The scheduler exposes **RPC methods** to:

  - Update runtime settings (e.g., ``max_calcjobs``, ``max_processes``).
  - Manipulate processes (e.g., ``play``, ``set_priority``).

Scheduler Management
--------------------

- Similar to the AiiDA daemon, the ``workgraph`` scheduler is managed using **``circus``**, providing reliable background process supervision and logging.
- Users can control and inspect active schedulers using the ``workgraph scheduler`` CLI:

  .. code-block:: bash

     workgraph scheduler start <name> [--max-calcjobs N] [--max-processes M]
     workgraph scheduler status <name>
     workgraph scheduler show <name>
     workgraph scheduler stop <name>
     workgraph scheduler list

Example
=======

1. Start the Scheduler
----------------------

Start a ``workgraph`` scheduler with limits on the number of concurrently running ``CalcJobs`` and total processes:

.. code-block:: bash

   workgraph scheduler start test-scheduler --max-calcjobs 2 --max-processes 10

2. Verify Scheduler Status
--------------------------

Check that the scheduler is running:

.. code-block:: bash

   workgraph scheduler status test-scheduler

Sample output:

.. code-block:: text

   Name            status    pk     waiting  running  calcjob  max_calcjobs  max_processes
   test-scheduler  Running   72507       0        0        0            2           10

3. Submit WorkGraphs with ArithmeticAddCalculation
--------------------------------------------------

Submit multiple ``WorkGraph`` instances, each containing several ``ArithmeticAddCalculation`` jobs. Each ``CalcJob`` is configured to sleep for 10 seconds to simulate runtime and help visualize scheduling limits.

.. code-block:: python

   from aiida_workgraph import WorkGraph
   from aiida import load_profile, orm
   from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

   load_profile()
   code = orm.load_code("add@localhost")  # Ensure this code is available

   for _ in range(4):  # Submit 4 WorkGraphs
       wg = WorkGraph("test_max_number_jobs")
       for i in range(5):  # Each WorkGraph has 5 CalcJobs
           task = wg.add_task(
               ArithmeticAddCalculation,
               name=f"add{i}",
               x=1,
               y=1,
               code=code
           )
           task.set({"metadata.options.sleep": 10})  # Simulate job runtime
       wg.submit(scheduler="test-scheduler")

4. Monitor Scheduler Progress
-----------------------------

Use the ``show`` command to inspect scheduler activity and confirm that job concurrency respects the specified limits:

.. code-block:: bash

   workgraph scheduler show test-scheduler

Sample output:

.. code-block:: text

   Report: Scheduler: test-scheduler
      PK    Created    Process label                    Process State   Priorities
   ------- ----------  -------------------------------  -------------- ------------
   72508   14s ago     WorkGraph<test_max_number_jobs>  ⏵ Waiting
   72509   13s ago     WorkGraph<test_max_number_jobs>  ⏵ Waiting
   72510   11s ago     WorkGraph<test_max_number_jobs>  ⏹ Created             -2
   ...
   72535   6s ago      ArithmeticAddCalculation         ⏹ Created             -1
   72538   5s ago      ArithmeticAddCalculation         ⏹ Created             -1

   Total results: 14

   name: test-scheduler
   pk: 72507
   running_process: 4
   waiting_process: 10
   running_calcjob: 2
   max_calcjobs: 2
   max_processes: 10

You can see that only 2 ``CalcJobs`` are running at a time (as per ``max_calcjobs=2``), and no more than 10 processes are handled concurrently.
