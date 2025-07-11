"""
========================
Grouping Tasks with Zone
========================
"""

# %%
## Introduction
# =============
# Zone task groups a collection of tasks. This will affect the execution order.
# - The zone and all its tasks are ready to run only if all links passing into the Zone are ready.
# - Tasks wait for the task inside the Zone needs to wait for the whole Zone to be finished.

#%%
# Setting up AiiDA environment
# ----------------------------
from aiida_workgraph import WorkGraph, task
from aiida import load_profile
load_profile()

#%%
# How to group tasks with `Zone`
# Here is an example of a Zone task:
# """


@task.calcfunction()
def add(x, y): return x + y

wg = WorkGraph("test_zone")
wg.ctx = {}
add1 = wg.add_task(add, name="add1", x=1, y=1)
zone1 = wg.add_task("workgraph.zone", name="zone1")
zone1.add_task(add, name="add2", x=1, y=1)
add3 = zone1.add_task(add, name="add3", x=1, y=add1.outputs.result)
zone1.add_task(add, name="add4", x=1, y=add3.outputs.result)
wg.add_task(add, name="add5", x=1, y=add3.outputs.result)
# export the workgraph to html file so that it can be visualized in a browser
wg.to_html()
# comment out the following line to visualize the workgraph in jupyter-notebook
# wg

###############################################################################
# ### Submit the WorkGraph and print out the report
#

from aiida.cmdline.utils.common import get_workchain_report
wg.submit(wait=True)
report = get_workchain_report(wg.process, "REPORT")

for line in report.split("\n"):
    if "tasks ready to run" in line:
        print(line)

###############################################################################
# We can see that the Zone is a collection of tasks. The Zone is ready to run only if all links passing into the Zone are ready. Tasks need input from the task inside the Zone need to wait for the whole Zone to be finished.
