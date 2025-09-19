"""
Run calculations remotely
=========================

"""

# %%
# .. _remote_calculations:
# Introduction
# ------------
#
# There are three ways to run calculations remotely using AiiDA:
#
# - `CalcJob`: run a shell command (via dedicated AiiDA plugins to manage inputs, outputs) on a remote computer (see :doc:`/howto/autogen/interoperate_with_aiida_core`).
# - `ShellJob`: run any shell command on a remote computer (see :doc:`/howto/autogen/shelljob`).
# - `PythonJob`: run any Python function on a remote computer.
#
# In this tutorial, we introduce how to run `PythonJob` tasks inside a WorkGraph.
# For more details on the `aiida-pythonjob` package itself, please refer to the `AiiDA-PythonJob documentation <https://aiida-pythonjob.readthedocs.io/en/latest/autogen/pythonjob.html>`_.

from aiida import load_profile

from aiida_workgraph import task

load_profile()

# %%
# PythonJob
# ---------
#
# One can create a `PythonJob` task using the ``@task.pythonjob`` decorator.
# This allows you to define a Python function that can be executed as a job on a remote computer.


@task.pythonjob
def add(x: int, y: int) -> int:
    return x + y


@task.graph
def RemoteAdd(x: int, y: int, computer: str) -> int:
    return add(x=x, y=y, computer=computer).result


wg = RemoteAdd.build(x=1, y=2, computer='localhost')
wg.run()

print('\nResult: ', wg.outputs.result.value)

# %%
# Understanding inputs and outputs
# --------------------------------
#
# When you apply the ``@task.pythonjob`` decorator, your Python function (like ``add``) transforms into an AiiDA `PythonJob` task.
# This task requires and produces more than just your function's direct inputs and outputs:
#
# * **Inputs:** In addition to your function's arguments (e.g., ``x``, ``y``), the task takes additional inputs (e.g., ``computer``) to manage its execution on a remote machine.
#
# * **Outputs:** In addition to your function's return value (available as ``.result`` by default), the task provides comprehensive outputs, such as ``remote_folder`` (a reference to the job's directory on the remote computer).

# %%
# Prepare Python environment on remote computer
# ---------------------------------------------
#
# `PythonJob` requires that the Python version on the remote computer matches the one used on the localhost where AiiDA is running.
# You can use `conda` to create a virtual environment with the same Python version, then activate the environment in the `metadata` of the `PythonJob` task.
#
# Here's an example demonstrating how to submit a `PythonJob` and chain it with another Python function within a `WorkGraph`.
# We'll also show how to pass custom scheduler commands via metadata, though for this example, they will be empty.


metadata = {
    'options': {
        # "custom_scheduler_commands" : "module load anaconda\nconda activate py3.11\n",
        'custom_scheduler_commands': '',  # Keeping it empty for this example
    }
}

# %%
# .. note::
#
#    Uncomment the ``custom_scheduler_commands`` line to modify it for your environment.

from typing import Annotated


@task
def multiply(x, y):
    return x * y


@task.graph
def RemoteAddLocalMultiply(x: int, y: int, computer: str, metadata: Annotated[dict, add.inputs.metadata]) -> int:
    the_sum = add(x=x, y=y, computer=computer, metadata=metadata).result
    return multiply(x=the_sum, y=4)  # this will run locally


wg = RemoteAddLocalMultiply.build(
    x=2,
    y=3,
    computer='localhost',
    metadata=metadata,
)
wg.run()

print('\nResult:', wg.outputs.result.value)


# %%
# Conclusion
# ----------
#
# In this tutorial, we briefly discussed the three ways to run calculations remotely using AiiDA: `CalcJob`, `ShellJob`, and `PythonJob`.
# We demonstrated how to use `PythonJob` to define and execute Python functions as remote jobs, including managing the remote Python environment via custom scheduler commands.
