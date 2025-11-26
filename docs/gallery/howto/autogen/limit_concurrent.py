"""
Limiting maximum running jobs in a WorkGraph
==============================================
"""

# %%
# Introduction
# ------------
#
# Imagine a workflow that needs to launch dozens of calculations. While this parallelism is powerful on a cluster, running them all at once on a local machine can overwhelm your system.
#
# `aiida-workgraph` provides a simple way to control the number of concurrent jobs for a single workflow.
# This is controlled by a single attribute:
#
# - ``WorkGraph.max_number_jobs``: Sets the maximum number of child processes that this *specific* WorkGraph can run simultaneously.

# sphinx_gallery_start_ignore
from aiida_workgraph.utils.logging import set_aiida_loglevel

set_aiida_loglevel('REPORT')
# sphinx_gallery_end_ignore

from aiida import load_profile
from aiida_workgraph import task

load_profile()

# %%
# Example
# --------
#
# Let's demonstrate this by wrapping an AiiDA `CalcJob` (`ArithmeticAddCalculation`) as a task and submitting several instances.
# We'll set a concurrency cap to prevent more than two jobs from running at the same time.

from aiida.orm import Int, load_code
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

# First, wrap the CalcJob class into a reusable WorkGraph task
AddTask = task(ArithmeticAddCalculation)

code = load_code('add@localhost')


# Define the WorkGraph that launches multiple independent jobs
@task.graph
def many_adds(n: int, code):
    """Launches 'n' independent addition jobs."""
    for k in range(n):
        AddTask(
            x=Int(k),
            y=Int(1),
            code=code,
        )


# Build the graph for 5 jobs
wg = many_adds.build(n=5, code=code)

# Set the maximum number of concurrent jobs to 2
wg.max_number_jobs = 2

# Run the workflow
wg.run()

# %%
# When you execute this, you will see messages like:
#
# ``The maximum number of subprocesses has been reached: 2. Cannot launch the job: ArithmeticAddCalculation2.``
#
# Once one of the running jobs completes, another will start, maintaining a maximum of two concurrent jobs until all five have finished.
#
# .. tip::
#
#     * **Workflow-Specific Limit:** The `max_number_jobs` attribute only governs child processes created by *this specific* WorkGraph instance. It does not limit other AiiDA jobs or workflows you might be running.
#     * **Nested graph task:** To limit the number of concurrent jobs in a nested graph task, you can set the `max_number_jobs`  in the `@task.graph` decorator, e.g. `@task.graph(max_number_jobs = 2)`.
#
#     * **Concurrency vs. Order:** This controls **how many** jobs run, not **which ones** run first. To enforce a specific execution sequence (e.g., Task B must run after Task A), you must define dependencies between them, see :ref:`Control task execution order <task_execution_order>`
#
#     * **Use with HPC Schedulers:** On a cluster with a scheduler like Slurm or PBS, it's often best to let the cluster's software manage global resources. However, `max_number_jobs` is still very useful for local development.
#     * **AiiDA-Scheduler** (experimental): For powerful, system-wide control over your AiiDA daemon, the `AiiDA-Scheduler <https://github.com/aiidateam/aiida-scheduler>`_ plugin is the recommended tool.
#


# sphinx_gallery_start_ignore
set_aiida_loglevel('ERROR')
# sphinx_gallery_end_ignore
