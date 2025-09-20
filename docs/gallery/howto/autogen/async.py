"""
Run ``async`` functions as tasks
================================

"""

# %%
# Introduction
# ------------
#
# The ``task`` decorator allows for the integration of ``asyncio`` within tasks, letting users control asynchronous functions.

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel('REPORT')
# sphinx_gallery_end_ignore

import asyncio

from aiida import load_profile
from aiida_workgraph import task

load_profile()


# %%
# Use the ``async`` function to make it non‑blocking:


@task
async def add_async(x, y, time=5):
    print(f'Sleeping for {time} seconds...')
    await asyncio.sleep(time)
    return x + y


@task
def multiply(x, y):
    return x * y


@task.graph
def AwaitableSum(x, y):
    sum1 = add_async(x, y, time=5).result
    sum2 = add_async(x, y, time=5).result
    return multiply(sum1, sum2).result


AwaitableSum.run(1, 2)

# %%
# Note the order of execution.
# The ``async`` tasks run concurrently (both print the sleep message immediately).
# Even though each sleeps for 5 seconds, they both complete around the same time (note the timestamps).
# Since the ``multiply`` task depends on both, it waits for both to finish before executing (note the timestamps).


# %%
# Summary
# -------
#
# In this section, we've explored the ``task`` decorator for integrating asynchronous functions within tasks.

# sphinx_gallery_start_ignore
set_aiida_loglevel('ERROR')
# sphinx_gallery_end_ignore
