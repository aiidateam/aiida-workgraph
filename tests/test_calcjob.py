import pytest
from aiida_workgraph import WorkGraph, task


def test_create_task_from_calcJob(add_code) -> None:
    """Test creating a task from a CalcJob."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    AddTask = task()(ArithmeticAddCalculation)
    with WorkGraph() as wg:
        outputs1 = AddTask(x=2, y=3, code=add_code)
        outputs2 = AddTask(x=outputs1.sum, y=3, code=add_code)
        wg.run()
    assert outputs2.sum.value == 8


@pytest.mark.usefixtures("started_daemon_client")
def test_submit(wg_calcjob: WorkGraph) -> None:
    """Submit simple calcjob."""
    wg = wg_calcjob
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True)
    assert wg.tasks.add2.outputs.sum.value == 9
