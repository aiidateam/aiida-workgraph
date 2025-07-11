"""
Write error-resistant workflows
===============================
"""


######################################################################
# Introduction
# ------------
# 
# In this tutorial, we will show how to implement the error handling in a
# WorkGraph. We will show how to implement the error handlers. They allow
# to execute a specific functions on a specific exit code of a
# ``CalcJob``. Every ``CalcJob`` defines its exit codes, on these exit
# codes we can define a function that is executed to handle the error and
# possibly restart the calculation.
# 
# Note that we error handling only works using by using concepts of the
# node graph programming. You might to first a look at the `node graph
# programming section <../node_graph_programming>`__.
# 


######################################################################
# Load the AiiDA environment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from aiida_workgraph import WorkGraph, task, Task

from aiida import load_profile, orm
from aiida.common.exceptions import NotExistent

load_profile()


try:
    bash_code = orm.load_code(
        "bash@localhost"
    )  # The computer label can also be omitted here
except NotExistent:
    bash_code = orm.InstalledCode(
        label="bash",
        computer=orm.load_computer("localhost"),
        filepath_executable="/bin/bash",
        default_calc_job_plugin="core.arithmetic.add",
    ).store()


######################################################################
# Exit codes
# ----------
# 
# We will take as example the ``ArithmeticAddCalculation`` from
# ``aiida-core`` CalcJob the methods we learn in this section can be
# applied to any ``CalcJob`` including also ``ShellJob`` and
# ``PythonJob``. If you are not familiar with the concept you might want
# to familiarize yourself in the `aiida-core
# documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/processes/concepts.html#topics-processes-concepts-exit-codes>`__
# Every ``CalcJob`` defines a serias of exit codes. we can print the exit
# codes simply by
# 

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
ArithmeticAddCalculation.exit_codes


######################################################################
# For this how-to we will write an error handler for exit code 410
# ``ERROR_NEGATIVE_NUMBER``.
# 

ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER


######################################################################
# Let us run a calculation that invokes this error code. If the computed
# sum of the inputs x and y is negative, the ``ArithmeticAddCalculation``
# fails with exit code 410.
# 

with WorkGraph("error_negative_number") as wg:
    wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code= bash_code)


wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)


######################################################################
# We can confirm that the task fails with this exit code in the CLI.
# 

# %verdi process status {wg.pk}


######################################################################
# Error handling
# --------------
# 
# To “register” a error handler for a WorkGraph, you simply define a
# function that takes the ``self`` and ``task_name`` as its arguments, and
# attach it as the ``error_hanlders`` of the WorkGraph.
# 
# You can specify the tasks and their exit codes that should trigger the
# error handler, as well as the maximum number of retries for a task:
# 
# .. code:: python
# 
#    tasks={"add1": {"exit_codes": [410],
#                    "max_retries": 5}
#          }
# 

def handle_negative_sum(task: Task):
    """Handle the failure code 410 of the `ArithmeticAddCalculation`.
    Simply make the inputs positive by taking the absolute value.
    """
    # modify task inputs
    task.set({"x": abs(task.inputs.x.value),
              "y": abs(task.inputs.y.value)})

with WorkGraph("handling_error_negative_number") as wg:
    wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code=bash_code)
    # Adding error handler logic
    wg.add_error_handler(handle_negative_sum, name="handle_negative_sum",
                         tasks={"add": {"exit_codes": [410], "max_retries": 5}})

wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)

# %verdi process status {wg.pk}


######################################################################
# Parametrized error handlers
# ---------------------------
# 
# One can also pass custom parameters to the error handler. For example,
# instead of simply make the inputs positive by taking the absolute value,
# we add an increment to the inputs. And the ``increment`` is a custom
# parameter of the error handler, which the user can specify when
# attaching the error handler to the WorkGraph, or update it during the
# execution of the WorkGraph.
# 


def handle_negative_sum(task: Task, increment: int = 1):
    """Handle the failure code 410 of the `ArithmeticAddCalculation`.
    Simply add an increment to the inputs.
    """
    # modify task inputs
    task.set({"x": task.inputs.x.value + increment,
              "y": task.inputs.y.value + increment})

with WorkGraph("handling_error_negative_number") as wg:
    wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code=bash_code)
    # Adding error handler logic
    wg.add_error_handler(handle_negative_sum, name="handle_negative_sum",
                         tasks={"add": {"exit_codes": [410],
                                        "max_retries": 5, # Note that retrying 5 times results in executing 6 times
                                        "kwargs": {"increment": 1}}})

wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)

# %verdi process status {wg.pk}


######################################################################
# We can confirm that the task first fails again with a 410. Then the
# WorkGraph restarts the task with the new inputs, and it finishes
# successfully.
# 

wg.reset()
wg.error_handlers.handle_negative_sum.tasks.add.kwargs.increment = 3
wg.run()

wg.reset()
wg.error_handlers['handle_negative_sum']['tasks']['add']['kwargs']['increment'] = 3
wg.run()


######################################################################
# Since we increase the inputs by a ``increment``, so it takes three
# retries before it finished successfully
# 

# %verdi process status {wg.pk}


######################################################################
# In this case, it only needs one retry to finish successfully.
# 


######################################################################
# Note that ``PythonJob`` task allows the user to attach the error handler
# directly to the task. Please check out the `aiida-pythonjob
# docs <https://aiida-pythonjob.readthedocs.io/en/latest/index.html>`__
# 


######################################################################
# Summary
# -------
# 


######################################################################
# Here we have shown how to implement error handling in a WorkGraph.
# 
