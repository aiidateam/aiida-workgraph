import time
import pytest
from aiida_workgraph import WorkGraph
from aiida.cmdline.utils.common import get_workchain_report


@pytest.mark.usefixtures("started_daemon_client")
def test_run_order(wg_engine: WorkGraph) -> None:
    """Test the order.
    Tasks should run in parallel and only depend on the input tasks."""
    wg = wg_engine
    wg.submit(wait=True)
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add0,add1" in report


@pytest.mark.skip(reason="The test is not stable.")
def test_reset_node(wg_engine: WorkGraph) -> None:
    """Modify a node during the excution of a WorkGraph."""
    wg = wg_engine
    wg.name = "test_reset"
    wg.submit()
    time.sleep(15)
    wg.tasks["add3"].set({"y": aiida.orm.Int(10).store()})
    wg.save()
    wg.wait()
    wg.update()
    assert wg.tasks["add5"].node.outputs.sum == 21
    assert wg.process.base.extras.get("_workgraph_queue_index") == 1
    assert len(wg.process.base.extras.get("_workgraph_queue")) == 1


@pytest.mark.usefixtures("started_daemon_client")
def test_max_number_jobs(add_code) -> None:
    from aiida_workgraph import WorkGraph
    from aiida.orm import Int
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    wg = WorkGraph("test_max_number_jobs")
    N = 3
    # Create N nodes
    for i in range(N):
        wg.add_task(
            ArithmeticAddCalculation, name=f"add{i}", x=Int(1), y=Int(1), code=add_code
        )
    # Set the maximum number of running jobs inside the WorkGraph
    wg.max_number_jobs = 2
    wg.submit(wait=True, timeout=100)
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add2" in report
