import pytest
from aiida_workgraph import WorkGraph, task
from aiida import orm


def test_run(wg_calcfunction: WorkGraph) -> None:
    """Run simple calcfunction."""
    wg = wg_calcfunction
    wg.name = "test_run_calcfunction"
    wg.run()
    print("state: ", wg.state)
    # print("results: ", results[])
    assert wg.tasks["sumdiff2"].node.outputs.sum == 9
    assert wg.tasks["sumdiff2"].outputs["sum"].socket_value == 9


@pytest.mark.usefixtures("started_daemon_client")
def test_dynamic_inputs() -> None:
    """Test dynamic inputs.
    For dynamic inputs, we allow the user to define the inputs manually.
    """

    @task.calcfunction(inputs=[{"name": "x"}, {"name": "y"}])
    def add(**kwargs):
        return kwargs["x"] + kwargs["y"]

    wg = WorkGraph("test_dynamic_inputs")
    wg.add_task(add, name="add1", x=orm.Int(1), y=orm.Int(2))
    wg.run()
    assert wg.tasks["add1"].outputs["result"].socket_value == 3
