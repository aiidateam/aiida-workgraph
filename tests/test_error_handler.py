from aiida_workgraph import WorkGraph
from aiida import load_profile, orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

load_profile()

code = orm.load_code("add@localhost")


def test_error_handlers():
    """Test error handlers."""
    from aiida.cmdline.utils.common import get_workchain_report

    def handle_negative_sum(self):
        """Handle negative sum by resetting the task and changing the inputs.
        self is the WorkGraph instance, thus we can access the tasks and the context.
        """
        node = self.get_task_state_info("add1", "process")
        if node and node.exit_code and node.exit_code.status == 410:
            self.report("Run error handler: handle_negative_sum.")
            self.reset_task("add1")
            # modify task inputs
            task = self.ctx.tasks["add1"]
            # the actual input values are stored in the properties dictionary
            task["properties"]["x"]["value"] = orm.Int(
                abs(-task["properties"]["x"]["value"])
            )
            task["properties"]["y"]["value"] = orm.Int(
                abs(-task["properties"]["y"]["value"])
            )

    wg = WorkGraph("restart_graph")
    wg.tasks.new(ArithmeticAddCalculation, name="add1")
    wg.error_handlers = [handle_negative_sum]
    wg.submit(
        inputs={
            "add1": {"code": code, "x": orm.Int(1), "y": orm.Int(-2)},
        },
        wait=True,
    )
    report = get_workchain_report(wg.process, "REPORT")
    assert "Run error handler: handle_negative_sum." in report
    assert wg.tasks["add1"].outputs["sum"].value == 3
