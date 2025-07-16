"""
Run tasks asynchronously
========================
"""

# %%
# Introduction
# ------------
#
# The ``awaitable`` decorator allows for the integration of ``asyncio`` within tasks, letting users control asynchronous functions.

# %%
# We'll temporarily set the AiiDA log level to ``REPORT``, so that we can inspect the execution order of the workgraph.

from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel("REPORT")

# %%
import asyncio

from aiida import load_profile
from aiida_workgraph import WorkGraph, task

load_profile()


# %%
@task.awaitable
async def awaitable_task(x, y):
    await asyncio.sleep(0.5)
    return x + y


with WorkGraph("AwaitableGraph") as wg:
    wg.inputs = dict.fromkeys(["x", "y"])
    awaitable_task(x=wg.inputs.x, y=wg.inputs.y)
    wg.inputs.x + wg.inputs.y

wg.run(
    inputs={
        "graph_inputs": {
            "x": 1,
            "y": 2,
        }
    },
)

# %%
# Note the timestamps. The addition task runs while the awaitable task sleeps.
# As the above tasks are functional, they would block one another.
# The ``awaitable`` decorator allows them to run concurrently.


# %%
# Notes on ``asyncio`` integration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The ``awaitable`` task lets the workgraph enter a ``WAITING`` state, yielding control to the ``asyncio`` event loop.
# This enables other tasks to run concurrently, though long-running calculations may delay the execution of awaitable tasks.

# %%
# Summary
# -------
#
# In this section, we've explored the ``awaitable`` decorator for integrating asynchronous functions within tasks.

set_aiida_loglevel("ERROR")
