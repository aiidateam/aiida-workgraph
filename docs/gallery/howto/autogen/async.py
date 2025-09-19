"""
Run ``async`` functions as tasks
================================

.. warning::
   **This feature is experimental.** The API for ``@task.awaitable`` is subject to change in future releases. We welcome your feedback on its functionality.

"""

# %%
# Introduction
# ------------
#
# The ``awaitable`` decorator allows for the integration of ``asyncio`` within tasks, letting users control asynchronous functions.

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel('REPORT')
# sphinx_gallery_end_ignore

import asyncio

from aiida import load_profile
from aiida_workgraph import task

load_profile()


# %%
# Use the ``@task.awaitable`` decorator on an ``async`` function to make it non‑blocking:


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


@task.awaitable
async def awaitable_add(x, y):
    await asyncio.sleep(4)
    return x + y


@task.graph
def AwaitableSum(x, y):
    async_sum = awaitable_add(x, y).result
    sync_sum = add(x, y).result
    return multiply(async_sum, sync_sum).result


wg = AwaitableSum.build(1, 2)
wg.run()


# %%
# Note the order of execution. The addition task runs while the awaitable task sleeps.
# As the above tasks are functional, they would block one another.
# The ``awaitable`` decorator allows them to run concurrently (both "ready to run").
# The ``add`` task finishes while the ``awaitable_add`` task sleeps.
# Since the ``multiply`` task depends on both, it waits for both to finish before executing (note the timestamps).


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

# sphinx_gallery_start_ignore
set_aiida_loglevel('ERROR')
# sphinx_gallery_end_ignore
