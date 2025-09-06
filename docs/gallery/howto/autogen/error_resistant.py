"""
Write error-resistant workflows
===============================
"""


######################################################################
# Introduction
# ------------
#
# In this how-to, we will show how to implement custom error handler functions that respond to exit codes with ``WorkGraph``, allowing you to automatically recover from errors, or to gracefully exit.
# We will walk through how to create error handlers that execute specific tasks based on the exit codes of ``CalcJob`` calculations.
# Error handling functionality requires using *nodegraph programming*.
# If you are unfamiliar with this approach, we recommend reviewing the `node graph programming section <../../concept/autogen/workgraph_concept.html>`__ first.
# We exemplify the error handling with the ``ArithmeticAddCalculation`` from ``aiida-core``, but the methods we learn in this section can be applied to any ``CalcJob``, including also ``ShellJob`` and ``PythonJob``.

######################################################################
# Load the AiiDA environment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from aiida_workgraph import WorkGraph, Task

from aiida import load_profile, orm
from aiida.common.exceptions import NotExistent
from aiida.cmdline.utils.ascii_vis import format_call_graph

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
# In AiiDA we can create error handlers for exit codes that are defined in the ``CalcJob`` specification.
# If you are not familiar with the concept of exit codes in ``CalcJobs`` s you might want to familiarize yourself in the `aiida-core documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/processes/concepts.html#topics-processes-concepts-exit-codes>`__
# For example for the ``ArithmeticAddCalculation`` we can the defined exit codes with
#

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

print(ArithmeticAddCalculation.exit_codes)


######################################################################
# We will write an error handler for the exit code 410 ``ERROR_NEGATIVE_NUMBER``.
#

print(ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER)


######################################################################
# If the computed sum of the inputs x and y is negative, the ``ArithmeticAddCalculation`` fails with exit code 410.
# We will run a calculation that exits with this error.
#

wg = WorkGraph("error_negative_number")
wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code=bash_code)


wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)


######################################################################
# We can confirm that the task fails with this exit code in the CLI.
#

print(format_call_graph(orm.load_node(wg.pk)))


######################################################################
# Error handling
# --------------
#
# To “register” a error handler for a WorkGraph, you simply define a function that takes the ``task`` as its arguments, and attach it as the ``error_handler`` of the ``WorkGraph``.
# You can specify the tasks and their exit codes that should trigger the error handler, as well as the maximum number of retries for a task:
#


def handle_negative_sum(task: Task):
    """Handle the failure code 410 of the `ArithmeticAddCalculation`.
    Simply make the inputs positive by taking the absolute value.
    """
    # modify task inputs
    task.set_inputs({"x": abs(task.inputs.x.value), "y": abs(task.inputs.y.value)})


wg = WorkGraph("handling_error_negative_number")
wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code=bash_code)
# Adding error handler logic
wg.add_error_handler(
    handle_negative_sum,
    name="handle_negative_sum",
    tasks={"add": {"exit_codes": [410], "max_retries": 5}},
)

wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)

######################################################################
#

print(format_call_graph(orm.load_node(wg.pk)))


######################################################################
# Parametrized error handlers
# ---------------------------
#
# We can also pass custom parameters to the error handler.
# For example, instead of simply making the inputs positive by taking the absolute value, we increment by an integer.
# The integer ``increment`` is a custom parameter of the error handler, which the user can specify when attaching the error handler to the WorkGraph.
#


def handle_negative_sum(task: Task, increment: int = 1):
    """Handle the failure code 410 of the `ArithmeticAddCalculation`.
    Simply add an increment to the inputs.
    """
    # modify task inputs
    task.set_inputs(
        {"x": task.inputs.x.value + increment, "y": task.inputs.y.value + increment}
    )


wg = WorkGraph("handling_error_negative_number")
wg.add_task(ArithmeticAddCalculation, name="add", x=1, y=-6, code=bash_code)
# Adding error handler logic
wg.add_error_handler(
    handle_negative_sum,
    name="handle_negative_sum",
    tasks={
        "add": {
            "exit_codes": [410],
            "max_retries": 5,  # Note that retrying 5 times results in executing 6 times
            "kwargs": {"increment": 1},
        }
    },
)

wg.run()
print("Task finished OK?:", wg.tasks.add.process.is_finished_ok)
print("Exit code        :", wg.tasks.add.process.exit_code)
print("Exit Message:    :", wg.tasks.add.process.exit_message)

######################################################################
#

print(format_call_graph(orm.load_node(wg.pk)))


######################################################################
# We can confirm that the task first fails again with a 410.
# Then the WorkGraph restarts the task with the new inputs, and it finishes successfully.
#

######################################################################
# We can also update the arguments of the error handler.
# Let us update the ``increment`` argument to 3 and restart the ``WorkGraph``.
#

# reset workgraph to start from the beginning
wg.reset()
wg.error_handlers["handle_negative_sum"]["tasks"]["add"]["kwargs"]["increment"] = 3
wg.run()


######################################################################
# In this case, it only needs one retry to finish successfully as adding two times 3 makes the ``y`` parameter positive
#

print(format_call_graph(orm.load_node(wg.pk)))


######################################################################
# .. note::
#
#    ``PythonJob`` task allows the user to attach the error handler directly to the task.
#    Please check out the `aiida-pythonjob documentation <https://aiida-pythonjob.readthedocs.io/en/stable/index.html>`__
#
