import pytest
from aiida_workgraph import WorkGraph, Task
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


@pytest.mark.usefixtures("started_daemon_client")
def test_error_handlers(add_code):
    """Test error handlers."""
    from aiida.cmdline.utils.common import get_workchain_report

    def handle_negative_sum(task: Task):
        """Handle negative sum by resetting the task and changing the inputs.
        self is the WorkGraph instance, thus we can access the tasks and the context.
        """
        # modify task inputs
        task.set(
            {
                "x": orm.Int(abs(task.inputs["x"].value)),
                "y": orm.Int(abs(task.inputs["y"].value)),
            }
        )
        msg = "Run error handler: handle_negative_sum."
        return msg

    wg = WorkGraph("restart_graph")
    wg.add_task(ArithmeticAddCalculation, name="add1")
    wg.add_error_handler(
        handle_negative_sum,
        name="handle_negative_sum",
        tasks={"add1": {"exit_codes": [410], "max_retries": 5, "kwargs": {}}},
    )
    wg.submit(
        inputs={
            "add1": {"code": add_code, "x": orm.Int(1), "y": orm.Int(-2)},
        },
        wait=True,
    )
    report = get_workchain_report(wg.process, "REPORT")
    assert "Run error handler: handle_negative_sum." in report
    assert wg.tasks["add1"].outputs["sum"].value == 3
