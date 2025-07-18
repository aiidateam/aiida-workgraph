"""

===============================
Run calculations remotely
===============================

"""

# %%
# Introduction
# ============
# There are three ways to run calculations remotely using AiiDA:
#
# - `CalcJob`: run a shell command (via dedicated AiiDA plugins to manage inputs, outputs) on a remote computer, please refer to `How to use aiida-core components inside WorkGraph <../../howto/autogen/use_calcjob_workchain.html>`_.
# - `ShellJob`: run any shell command on a remote computer, please refer to the `Run shell commands as a task <../../howto/autogen/shelljob.html>`_.
# - `PythonJob`: run any Python function on a remote computer.
#
# In this tutorial, we introduce how to run `PythonJob` tasks inside a WorkGraph.
# For more details on the `aiida-pythonjob` package itself, please refer to the `AiiDA-PythonJob documentation <https://aiida-pythonjob.readthedocs.io/en/latest/autogen/pythonjob.html>`_.

from aiida import load_profile

load_profile()

# %%
# PythonJob task
# ========================
# One can create a `PythonJob` task using the ``@task.pythonjob()`` decorator.
# This allows you to define a Python function that can be executed as a job on a remote computer.

from aiida_workgraph import WorkGraph, task

# Decorator to define a pythonjob
@task.pythonjob()
def add(x, y):
    return x + y


with WorkGraph() as wg:
    # Call the add task, specifying 'localhost' as the computer.
    # Notice we pass 'computer="localhost"' even though it's not
    # an argument to the 'add(x, y)' function directly.
    outputs = add(x=1, y=2, computer="localhost")
    wg.run()

# Print out the result once the WorkGraph has finished executing.
print("\nResult: ", outputs.result)

# %%
# Understanding Inputs and Outputs
# --------------------------------
# When you apply the ``@task.pythonjob()`` decorator, your Python function (like ``add``) transforms into an AiiDA `PythonJob` task.
# This task requires and produces more than just your function's direct inputs and outputs:
#
# * **Inputs:** Besides your function's arguments (e.g., ``x``, ``y``), the task takes additional inputs like ``computer`` to manage its execution on a remote machine.
#
# * **Outputs:** Beyond your function's return value (available as ``.result`` by default), the task provides comprehensive outputs, such as ``remote_folder`` (a reference to the job's directory on the remote computer).

# %%
# Prepare Python environment on remote computer
# -----------------------------------------------
# `PythonJob` requires that the Python version on the remote computer matches the one used on the localhost where AiiDA is running.
# You can use `conda` to create a virtual environment with the same Python version, then activate the environment in the `metadata` of the `PythonJob` task.
#
# Here's an example demonstrating how to submit a `PythonJob` and chain it with another Python function within a `WorkGraph`.
# We'll also show how to pass custom scheduler commands via metadata, though for this example, they will be empty.
#


# Uncomment the 'custom_scheduler_commands' line below and modify it to suit your environment.
metadata = {
    "options": {
        # 'custom_scheduler_commands' : 'module load anaconda\nconda activate py3.11\n',
        "custom_scheduler_commands": "",  # Keeping it empty for this example
    }
}

# This is a normal task. It will be executed locally within the WorkGraph.
@task()
def multiply(x, y):
    return x * y


with WorkGraph(name="test_pythonjob") as wg:
    # We pass the metadata to ensure any custom scheduler commands are applied.
    outputs1 = add(x=2, y=3, computer="localhost", metadata=metadata)
    # Add another task which will run locally
    wg.outputs.result = multiply(x=outputs1.result, y=4).result
    wg.run()

# ------------------------- Print the output -------------------------
print("\nResult {} \n\n".format(wg.outputs.result.value))


# %%
# Conclusion
# ===========
# In this tutorial, we demonstrated the three ways to run calculations remotely using AiiDA:
# `CalcJob`, `ShellJob`, and `PythonJob`. We focused on `PythonJob`, showing how to
# define and execute Python functions as remote jobs,
# including managing the remote Python environment via custom scheduler commands.
