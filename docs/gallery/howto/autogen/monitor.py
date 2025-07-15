"""
Assign a monitoring task
========================
"""


# %%
# Introduction
# ------------
#
# The ``monitor`` decorator is designed for tasks that need to poll a specific state at regular intervals until a success criterion is met. This is useful for various scenarios, including time-based triggers, file existence checks, and monitoring other tasks or workgraphs.
#
# Possible use Cases
# ~~~~~~~~~~~~~~~~~~
#
# - **Time-based events**: Start a task at a specified time
# - **File-based events**: Execute a task when a particular file exists
# - **Monitor a task**: Observe the state of another task and act based on certain conditions
# - **Cross-workGraph dependencies**: Check the state of a task in a different workgraph
#
# Behavior
# ~~~~~~~~
#
# While polling, the task sleeps for a specified interval (default 1.0 second), allowing the workgraph engine to manage other tasks.
#
# Example usage
# ~~~~~~~~~~~~~
#
# The monitor task has two built-in parameters:
#
# - ``interval``: The time interval between each poll
# - ``timeout``: The maximum time to wait for the success criterion to be met
#
# In the following sections, we will walk through examples of how to use the `monitor` task decorator for these scenarios.
#

from aiida_workgraph import WorkGraph, task
from aiida import load_profile

_ = load_profile()


# %%
# Time-based events
# -----------------
#
# Here we design a workgraph that waits until 5 seconds have passed.

import datetime


@task.monitor
def time_monitor(time):
    return datetime.datetime.now() > datetime.datetime.fromisoformat(time.value)


with WorkGraph("TimeMonitor") as wg:
    wg.inputs = {
        "monitor": {"time": None},
        "add": dict.fromkeys(["x", "y"]),
    }
    time_monitor(time=wg.inputs.monitor.time) >> wg.inputs.add.x + wg.inputs.add.y

wg

wg.run(
    inputs={
        "graph_inputs": {
            "monitor": {
                "time": (
                    datetime.datetime.now() + datetime.timedelta(seconds=5)
                ).isoformat(),
            },
            "add": {
                "x": 1,
                "y": 2,
            },
        }
    }
)

# Note the time difference between the monitor task and the next (~5 seconds)


# %%
# File-based events
# -----------------
#
# Here we design a workgraph that waits until a file exists before it proceeds. We create the file asynchronously after 5 seconds. During this time, the monitor task will poll for the file's existence. Once the file is found, the workgraph will proceed to discard the file and add two numbers (independently), then finally multiply the sum by a factor.

import os
import asyncio
from node_graph.group import group


@task.awaitable
async def sleep_create_file(filepath, content):
    await asyncio.sleep(5)
    with open(filepath, "w") as f:
        f.write(content)


@task.monitor
def file_monitor(filepath):
    return os.path.exists(filepath)


@task
def discard_file(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
    else:
        raise FileNotFoundError(f"File {filepath} does not exist.")


with WorkGraph("FileMonitor") as wg:
    wg.inputs = {
        "add": dict.fromkeys(["x", "y"]),
        "multiply": dict.fromkeys(["factor"]),
    }

    # Asynchronously create a file after 5 seconds
    sleep_create_file(
        filepath="/tmp/monitor_test.txt",
        content="This is a test file",
    )

    # While above is running, monitor for the file
    # Once done (`>>` means wait), add two numbers and discard the file
    # Finally, multiply the sum by a factor
    (
        file_monitor(
            filepath="/tmp/monitor_test.txt",
            interval=1,
            timeout=10,
        )
        >> group(
            discard_file(filepath="/tmp/monitor_test.txt"),
            the_sum := wg.inputs.add.x + wg.inputs.add.y,
        )
        >> the_sum * wg.inputs.multiply.factor
    )

wg

wg.run(
    inputs={
        "graph_inputs": {
            "add": {
                "x": 1,
                "y": 2,
            },
            "multiply": {
                "factor": 3,
            },
        }
    }
)

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
#    monitor1 = wg.add_task("workgraph.time_monitor", datetime=datetime.datetime.now() + datetime.timedelta(seconds=10))
#    monitor2 = wg.add_task("workgraph.file_monitor", filepath="/tmp/test.txt")

# %%
# Summary
# -------
#
# You have learned how to use the ``monitor`` decorator to create tasks that poll for specific conditions, such as time-based events, file-based events, and task monitoring. You also learned how to kill a monitor task and about the built-in monitor tasks provided by `WorkGraph`.
