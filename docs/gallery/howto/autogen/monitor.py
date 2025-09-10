"""
Monitor external events as a task
=================================
"""


# %%
# Introduction
# ------------
#
# The ``monitor`` decorator is designed for tasks that need to poll a specific state at regular intervals until a success criterion is met.
# Some possible use cases include:
#
# - **Time-based events**: Start a task at a specified time
# - **File-based events**: Execute a task when a particular file exists
# - **Monitor a task**: Observe the state of another task and act based on certain conditions
# - **Cross-workGraph dependencies**: Check the state of a task in a different workgraph
#
# While polling, the task sleeps for a specified interval (default 1.0 second), allowing the workgraph engine to manage other tasks.
#
# The monitor task has two built-in parameters:
#
# - ``interval``: The time interval between each poll
# - ``timeout``: The maximum time to wait for the success criterion to be met
#
# In the following sections, we will walk through examples of how to use the `monitor` task decorator for these scenarios.

# %%
# We'll temporarily set the AiiDA log level to ``REPORT``, so that we can inspect the execution order of the workgraph.

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel("REPORT")
# sphinx_gallery_end_ignore

from aiida_workgraph import task
from aiida import load_profile

load_profile()

# %%
# Time-based events
# -----------------
#
# Here we design a workgraph that waits until 5 seconds have passed.

import datetime


@task.monitor
def monitor_time(time: datetime.datetime):
    return datetime.datetime.now() > time.value


@task
def add(x, y):
    return x + y


@task.graph
def TimeMonitor(time, x, y):
    monitor_time(time) >> (the_sum := add(x, y).result)  # wait, THEN (>>) add
    return the_sum


wg = TimeMonitor.build_graph(
    time=datetime.datetime.now() + datetime.timedelta(seconds=5),
    x=1,
    y=2,
)
wg.to_html()

# %%
wg.run()

# %%
# Note the time difference between the monitor task and the next (~5 seconds)

# %%
# File-based events
# -----------------
#
# Here we design a workgraph that waits until a file exists before it proceeds. We create the file asynchronously after 5 seconds. During this time, the monitor task will poll for the file's existence. Once the file is found, the workgraph will proceed to discard the file and add two numbers (independently), then finally multiply the sum by a factor.

import os
import asyncio
from aiida_workgraph.collection import group


@task.awaitable
async def sleep_create_file(filepath, content):
    await asyncio.sleep(5)
    with open(filepath, "w") as f:
        f.write(content)


@task.monitor
def monitor_file(filepath):
    return os.path.exists(filepath)


@task
def discard_file(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
    else:
        raise FileNotFoundError(f"File {filepath} does not exist.")


@task
def multiply(x, y):
    return x * y


@task.graph
def FileMonitor(x, y, z):
    # Asynchronously create a file after 5 seconds
    sleep_create_file(
        filepath="/tmp/monitor_test.txt",
        content="This is a test file",
    )

    # While above is running, monitor for the file
    # Once done (`>>` means wait on left operant), discard the file and add two numbers
    monitor_file("/tmp/monitor_test.txt", interval=1, timeout=10,) >> group(
        discard_file("/tmp/monitor_test.txt"),
        the_sum := add(x, y).result,
    )

    # Finally, multiply the sum by a factor
    return multiply(the_sum, z).result


wg = FileMonitor.build_graph(x=1, y=2, z=3)
wg.to_html()

# %%
wg.run()

# %%
# You can inspect the process reports above to verify the order of events.

# %%
# Kill a monitor task
# -------------------
#
# One can kill a running monitor task by using the following command:
#
# .. code:: console
#
#    workgraph task kill <workgraph_pk> <task_name>
#
# For example:
#
# .. code:: console
#
#    workgraph task kill 119974 monitor1
#
# A killed task will has the status ``KILLED`` and the following task will not be executed.
#

# %%
# Built-in monitors
# -----------------
#
# ``WorkGraph`` provides the the above time and file monitors as built-in tasks. You can use them directly without having to define them.
#
# .. code:: python
#
#    from aiida_workgraph.tasks.monitors import monitor_time, monitor_file

# %%
# Summary
# -------
#
# You have learned how to use the ``monitor`` decorator to create tasks that poll for specific conditions, such as time-based events, file-based events, and task monitoring. You also learned how to kill a monitor task and about the built-in monitor tasks provided by `WorkGraph`.

# sphinx_gallery_start_ignore
set_aiida_loglevel("ERROR")
# sphinx_gallery_end_ignore
