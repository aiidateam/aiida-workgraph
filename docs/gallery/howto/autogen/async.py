"""
Run tasks asynchronously
========================
"""

# %%
# Introduction
# ------------
#
# The ``awaitable`` decorator allows for the integration of ``asyncio`` within tasks, letting users control asynchronous functions.

from aiida_workgraph import WorkGraph, task
from aiida import load_profile

_ = load_profile()

import asyncio


@task.awaitable
async def awaitable_task(x, y):
    await asyncio.sleep(0.5)
    return x + y


with WorkGraph("AwaitableGraph") as wg:
    wg.inputs = dict.fromkeys(["x", "y"])
    awaitable_task(x=wg.inputs.x, y=wg.inputs.y)

wg.run(
    inputs={
        "graph_inputs": {
            "x": 1,
            "y": 2,
        }
    }
)

# %%
# Notes on asyncio Integration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The awaitable task lets the workgraph enter a ``Waiting`` state, yielding control to the asyncio event loop. This enables other tasks to run concurrently, though long-running calculations may delay the execution of awaitable tasks.

# %%
# Summary
# -------
#
# In this section, we've explored the ``awaitable`` decorator for integrating asynchronous functions within tasks.
