from aiida_workgraph import WorkGraph
from aiida import load_profile, orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

load_profile()

code = orm.load_code("add@localhost")


def test_error_handlers():
    """Test error handlers."""
    from aiida.cmdline.utils.common import get_workchain_report

    def handle_negative_sum(self, task_name: str, **kwargs):
        """Handle negative sum by resetting the task and changing the inputs.
        self is the WorkGraph instance, thus we can access the tasks and the context.
        """
        self.report("Run error handler: handle_negative_sum.")
        task = self.get_task(task_name)
        # modify task inputs
        task.set(
            {
                "x": orm.Int(abs(task.inputs["x"].value)),
                "y": orm.Int(abs(task.inputs["y"].value)),
            }
        )
        self.update_task(task)

    wg = WorkGraph("restart_graph")
    wg.tasks.new(ArithmeticAddCalculation, name="add1")
    wg.attach_error_handler(
        handle_negative_sum,
        name="handle_negative_sum",
        tasks={"add1": {"exit_codes": [410], "max_retries": 5, "kwargs": {}}},
    )
    wg.submit(
        inputs={
            "add1": {"code": code, "x": orm.Int(1), "y": orm.Int(-2)},
        },
        wait=True,
    )
    report = get_workchain_report(wg.process, "REPORT")
    assert "Run error handler: handle_negative_sum." in report
    assert wg.tasks["add1"].outputs["sum"].value == 3
