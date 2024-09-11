"""
=====================
Run tasks in parallel
=====================
"""
# %%
# Introduction
# ============
# In this tutorial, you will learn how to run tasks and WorkGraphs in parallel.
# When defining the dependencies WorkGraph by linking tasks the WorkGraph
# engine will automatically take care of parallelizing the independent tasks. One
# caveat is that we cannot use calcfunctions for this purpose as they all run
# in the same runner environment and therefore are blocking each other. For
# that reason we need to use `CalcJob`s that can be run in different runner
# environments and therefore can be run in parallel.

# Load the AiiDA profile.
from aiida import load_profile

load_profile()


# %%
# Parallel addition workflow
# ==========================
# Suppose we want to calculate ```x + y + u + v``` in a parallel, instead of
# computing sequentially ```(((x + y) + u) + v)``` we compute it like
# ```((x + y) + (u + v))``` to compute ```x + y``` and ```u + v``` in parallel.
# aiida-core already provides a ArithmeticAddCalculation CalcJob for performing
# addition which we will use it for this example

from aiida_workgraph import WorkGraph, task
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida.orm import InstalledCode, load_computer, load_code, load_node
from aiida.common.exceptions import NotExistent

# The ArithmeticAddCalculation needs to know where bash is stored
try:
    code = load_code("add@localhost")  # The computer label can also be omitted here
except NotExistent:
    code = InstalledCode(
        computer=load_computer("localhost"),
        filepath_executable="/bin/bash",
        label="add",
        default_calc_job_plugin="core.arithmetic.add",
    ).store()

wg = WorkGraph("parallel")
x, y, u, v = (1, 2, 3, 4)
add_xy = wg.add_task(ArithmeticAddCalculation, name="add_xy", x=x, y=y, code=code)
add_xy.set({"metadata.options.sleep": 3})  # the CalcJob will sleep 3 seconds
add_uv = wg.add_task(ArithmeticAddCalculation, name="add_uv", x=u, y=v, code=code)
add_uv.set({"metadata.options.sleep": 3})  # the CalcJob will sleep 3 seconds
add_xyuv = wg.add_task(
    ArithmeticAddCalculation,
    name="add_xyuv",
    x=add_xy.outputs["sum"],
    y=add_uv.outputs["sum"],
    code=code,
)
# %%
# We can verify that the tasks add_xy and add_uv are independent from each other
# and therefore will be run automatically in parallel.

wg.to_html()


# %%
# Running workgraph

wg.submit(wait=True)

# %%
#  We look at the ctime (the time of creation when submitted/run) and the mtime (the time the task has been last modified which is when its state changes to finish).
print("add_xy created at:", add_xy.ctime.time(), "finished at:", add_xy.mtime.time())
print("add_uv created at:", add_uv.ctime.time(), "finished at:", add_uv.mtime.time())

# %%
# We can see that both CalcJob's have been created almost at the same time

# %%
# Comparison with a calcfunction
# ------------------------------
#


@task.calcfunction()
def add(x, y, sleep):
    import time

    time.sleep(sleep.value)
    return x + y


wg = WorkGraph("parallel")
x, y, u, v = (1, 2, 3, 4)
add_xy = wg.add_task(add, x=x, y=y, sleep=3)
add_uv = wg.add_task(add, x=x, y=y, sleep=3)
add_xyuv = wg.add_task(
    add, x=add_xy.outputs["result"], y=add_uv.outputs["result"], sleep=0
)

wg.to_html()

# %%

wg.submit(wait=True)

# %%
# Printing timings

print("add_xy created at", add_xy.ctime.time(), "finished at", add_xy.mtime.time())
print("add_uv created at", add_uv.ctime.time(), "finished at", add_uv.mtime.time())

# %%
# We can see that the calcfunctions have been run with a 3 seconds delay


# %%
# Parallelizing WorkGraphs
# ========================
# We will parallelize a workgraph by two ways, one time we submit all workgraphs,
# the other time we use the graph builder to submit once the whole workflow.


# This is our initial WorkGraph we want to parallelize
@task.graph_builder(
    inputs=[{"name": "integer"}], outputs=[{"name": "sum", "from": "sum_task.result"}]
)
def add10(integer):
    wg = WorkGraph()
    code = load_code("add@localhost")  # code needs to loaded in the graph builder
    add = wg.add_task(
        ArithmeticAddCalculation, name="sum_task", x=10, y=integer, code=code
    )
    add.set({"metadata.options.sleep": 3})
    return wg


# %%

wgs = []
for i in range(2):
    wg = WorkGraph(f"parallel_wg{i}")
    wg.add_task(add10, name=f"add10_{i}", integer=i)
    wgs.append(wg)

# We use wait=False so we can continue submitting
wgs[0].submit()  # do not wait (by default), so that we can continue to submit next WG.
wgs[1].submit(wait=True)
# we wait for all the WorkGraphs to finish
wgs[0].wait()

# %%
# We print the difference between the mtime (the time the WorkGraph has been last time changed) and the ctime (the time of creation). Since the WorkGraph's status is changed when finished, this give us a good estimate of the running time.
print(
    "WG0 created:",
    load_node(wgs[0].pk).ctime.time(),
    "finished:",
    load_node(wgs[0].pk).mtime.time(),
)
print("Time WG0", load_node(wgs[0].pk).mtime - load_node(wgs[0].pk).ctime)
print(
    "WG1 created:",
    load_node(wgs[1].pk).ctime.time(),
    "finished:",
    load_node(wgs[1].pk).mtime.time(),
)
print("Time WG1", load_node(wgs[1].pk).mtime - load_node(wgs[1].pk).ctime)


# %%
# Using graph builder
# -------------------


# This graph_builder runs the add10 over a loop and its
@task.graph_builder()
def parallel_add(nb_iterations):
    wg = WorkGraph()
    for i in range(nb_iterations):
        wg.add_task(add10, name=f"add10_{i}", integer=i)
    return wg


# Submitting a parallel that adds 10 two times to different numbers
wg = WorkGraph(f"parallel_graph_builder")
add_task = wg.add_task(parallel_add, name="parallel_add", nb_iterations=2)
wg.to_html()

# %%
wg.submit(wait=True)


# %%
# We look at the times of creation and last change
print("Time for running with graph builder", add_task.mtime - add_task.ctime)

# %%
# We can see that the time is less than 6 seconds which means that the two additions
# were performed in parallel

# %%
# Increasing number of daemon workers
# -----------------------------------
# Since each daemon worker can only manage one WorkGraph (handling the results)
# at a time, one can experience slow downs when running many jobs that can be
# run in parallel. The optimal number of workers depends highly on the jobs
# that are run.

from aiida.engine.daemon.client import get_daemon_client

client = get_daemon_client()

# %%
# We rerun the last graph builder for 3 iterations

wg = WorkGraph("wg_daemon_worker_1")
wg.add_task(parallel_add, name="parallel_add", nb_iterations=3)
wg.to_html()

# %%
wg.submit(wait=True)
print(
    f"Time for running with {client.get_numprocesses()['numprocesses']} worker",
    load_node(wg.pk).mtime - load_node(wg.pk).ctime,
)

# %%
# We increase the number of workers by one. One can also do this in the workgraph GUI.

client = get_daemon_client()
client.increase_workers(1)

# %%
# Now we submit again and the time have shortens a bit.

wg = WorkGraph("wg_daemon_worker_2")
wg.add_task(parallel_add, name="parallel_add", nb_iterations=3)
wg.to_html()

# %%

wg.submit(wait=True)
print(
    f"Time for running with {client.get_numprocesses()['numprocesses']} worker",
    load_node(wg.pk).mtime - load_node(wg.pk).ctime,
)

# %%
# Note that on readthedocs you will not see a big difference due to the hardware.
# With a limited number of CPU the workers cannot be parallelized

import multiprocessing

print("Number of CPUs", multiprocessing.cpu_count())

# %%
# Reset back to one worker
client.decrease_workers(1)

# %%
# Maximum number of active WorkGraphs
# -----------------------------------
# Be aware that for the moment AiiDA can only run 200 WorkGraphs at the same time.
# To increase that limit one can set this variable to a higher value.
#
#     .. code-block:: bash
#
#        verdi config set daemon.worker_process_slots 200
#        verdi daemon restart


# %%
# Further reading
# ---------------
# Now you learned how to run tasks in parallel you might want to know how to
# aggregate the results of all these parallel tasks (e.g. taking the mean of
# all computed values). For this you can further read
# :ref:`sphx_glr_howto_autogen_aggregate.py`.
