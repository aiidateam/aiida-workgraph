"""
===========================
Error Handling in WorkGraph
===========================

This tutorial demonstrates how to implement robust error handling in an AiiDA WorkGraph.
We will compare the traditional approach using `BaseRestartWorkChain` with the more
flexible error handling mechanism available in `aiida-workgraph`.
"""
# %%
# Introduction: Error Handling in AiiDA
# =====================================
# AiiDA provides a ``BaseRestartWorkChain`` class to create workflows that can
# automatically recover from known failure modes of calculations. This is
# typically done by defining specific handlers for error codes produced by a
# particular ``Process`` class.
#
# While powerful, this approach ties the error handling logic to a specific
# WorkChain implementation. In contrast, `aiida-workgraph` offers a more general
# mechanism where error handlers can be defined as independent functions and
# attached to any task within a WorkGraph.
#
# First, let's look at a `BaseRestartWorkChain` example.

from aiida.engine import BaseRestartWorkChain
from aiida.plugins import CalculationFactory
from aiida import orm
from aiida.engine import while_
from aiida.engine import process_handler, ProcessHandlerReport

ArithmeticAddCalculation = CalculationFactory('core.arithmetic.add')


class ArithmeticAddBaseWorkChain(BaseRestartWorkChain):
    _process_class = ArithmeticAddCalculation

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.expose_inputs(ArithmeticAddCalculation, namespace='add')
        spec.expose_outputs(ArithmeticAddCalculation)
        spec.outline(
            cls.setup,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the process in the
        internal loop.
        """
        super().setup()
        self.ctx.inputs = self.exposed_inputs(ArithmeticAddCalculation, 'add')

    @process_handler
    def handle_negative_sum(self, node):
        """Check if the calculation failed with `ERROR_NEGATIVE_NUMBER`.

        If this is the case, simply make the inputs positive by taking the absolute value.

        :param node: the node of the subprocess that was ran in the current iteration.
        :return: optional :class:`~aiida.engine.processes.workchains.utils.ProcessHandlerReport` instance to signal
            that a problem was detected and potentially handled.
        """
        if node.exit_status == ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER.status:
            self.ctx.inputs['x'] = orm.Int(abs(node.inputs.x.value))
            self.ctx.inputs['y'] = orm.Int(abs(node.inputs.y.value))
            return ProcessHandlerReport()


# %%
# Error Handling in a WorkGraph
# =============================
# Now, let's implement the same logic using a WorkGraph. The process involves
# three main steps:
#
# 1.  **Define an error handler function**: This function takes a `Task` object
#     as input and contains the logic to resolve the error.
# 2.  **Build a WorkGraph**: Create a WorkGraph and add the task that might fail.
# 3.  **Register the error handler**: Attach the handler function to the specific
#     task, specifying which exit codes should trigger it.

from aiida_workgraph import WorkGraph, Task
from aiida.cmdline.utils.ascii_vis import format_call_graph


def handle_negative_sum(task: Task):
    """
    Error handler for ArithmeticAddCalculation.

    This function is triggered when the calculation fails with exit code 410
    (ERROR_NEGATIVE_NUMBER). It modifies the task's inputs by taking the
    absolute value of `x` and `y`, allowing the task to be restarted.

    :param task: The Task object that failed.
    :return: A message indicating that the handler was executed.
    """
    from aiida.orm import Int

    print(f'Executing error handler for task: {task.name}')
    # Modify the inputs of the failed task for the next retry.
    task.set_inputs({'x': Int(abs(task.inputs['x'].value)), 'y': Int(abs(task.inputs['y'].value))})
    msg = 'Run error handler: handle_negative_sum.'
    return msg


# %%
# Building and Running the WorkGraph
# ==================================
# We create a WorkGraph, add the `ArithmeticAddCalculation` as a task, and then
# register our `handle_negative_sum` function as its error handler.

# 1. Create a new WorkGraph
wg = WorkGraph('error_handling_graph')

# 2. Add the calculation task
task1 = wg.add_task(ArithmeticAddCalculation, name='add_task')

# 3. Register the error handler for the task
task1.add_error_handler(
    {
        'handle_negative_sum': {
            'executor': handle_negative_sum,
            'exit_codes': [ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER.status],
            'max_retries': 3,
        }
    }
)

# ------------------------- Submit the calculation -------------------
from aiida import load_profile
from aiida.orm import Int, load_code

load_profile()

# Define the inputs that will cause the calculation to fail initially.
inputs = {
    'add_task': {
        'code': load_code('add@localhost'),
        'x': Int(1),
        'y': Int(-6),  # This will result in a negative sum, triggering the error.
    },
}

# Submit the WorkGraph and wait for it to complete.
# The error handler will be executed automatically by the engine.
print('Submitting WorkGraph that is expected to fail and be corrected...')
wg.run(inputs=inputs)


# %%
# Verifying the Result
# ====================
# After the WorkGraph finishes, we can inspect the final state of the task
# to confirm that the error was handled and the calculation eventually
# succeeded.

add_task = wg.tasks['add_task']
print(f'Task finished OK? {add_task.process.is_finished_ok}')
print(f'Final exit status: {add_task.process.exit_status}')
print(f"Final result: {add_task.outputs['sum'].value}")

# We can also visualize the process call graph to see the failure and restart.
# The graph will show the first `ArithmeticAddCalculation` failing (exit code 410)
# and a second, corrected one finishing successfully (exit code 0).
print('\n--- Process Call Graph ---')
print(format_call_graph(wg.process))


# %%
# Conclusion
# ==========
# The `aiida-workgraph` error handling mechanism provides a powerful and
# decoupled way to manage failures. By separating the error handling logic
# from the workflow definition, you can create reusable handlers that can be
# applied to different tasks across various WorkGraphs, leading to more
# modular and maintainable scientific workflows.
#
# For detailed documentation, please refer to the `Write error-resistant workflows <../../howto/autogen/error_resistant>`_).
