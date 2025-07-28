from aiida_workgraph import WorkGraph, Task
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


def test_error_handlers(add_code, capsys):
    """Test error handlers."""

    def handle_negative_sum(task: Task):
        """Handle negative sum by resetting the task and changing the inputs.
        self is the WorkGraph instance, thus we can access the tasks and the context.
        """
        # modify task inputs
        task.set(
            {
                "x": orm.Int(abs(task.inputs.x.value)),
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
        tasks={
            "add1": {
                "exit_codes": [
                    ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER
                ],
                "max_retries": 5,
                "kwargs": {},
            }
        },
    )
    assert len(wg.error_handlers) == 1
    wg.run(
        inputs={
            "add1": {"code": add_code, "x": orm.Int(1), "y": orm.Int(-2)},
        },
    )
    captured = capsys.readouterr()
    report = captured.out
    assert "Run error handler: handle_negative_sum." in report
    assert wg.tasks.add1.outputs.sum.value == 3
