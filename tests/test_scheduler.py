import pytest
from aiida_workgraph import WorkGraph
from aiida_workgraph.tools.pytest_fixtures.scheduler import DEFAULT_TEST_SCHEDULER_NAME
from aiida_workgraph.engine.scheduler import Scheduler
from aiida_workgraph.engine.scheduler.client import get_scheduler
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
import datetime


@pytest.mark.usefixtures("started_daemon_client")
@pytest.mark.usefixtures("started_scheduler_client")
def test_submit(add_code) -> None:
    """Submit simple calcjob."""
    wg = WorkGraph(name="test_debug_math")
    add1 = wg.add_task(ArithmeticAddCalculation, "add1", x=2, y=3, code=add_code)
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=4, y=3, code=add_code)
    add1.set({"metadata.options.sleep": 5})
    add2.set({"metadata.options.sleep": 5})
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True, timeout=30, scheduler=DEFAULT_TEST_SCHEDULER_NAME)
    wg.update()
    assert add2.process.mtime - add1.process.mtime < datetime.timedelta(seconds=5)

    #
    status = Scheduler.get_status(name=DEFAULT_TEST_SCHEDULER_NAME)
    assert status is not None
    assert "running_process" in status
    status = Scheduler.set_max_calcjobs(
        name=DEFAULT_TEST_SCHEDULER_NAME, max_calcjobs=1
    )
    scheduler = get_scheduler(name=DEFAULT_TEST_SCHEDULER_NAME)
    assert scheduler.max_calcjobs == 1
    #
    wg.reset()
    wg.submit(wait=True, timeout=30, scheduler=DEFAULT_TEST_SCHEDULER_NAME)
    wg.update()
    assert add2.process.mtime - add1.process.mtime > datetime.timedelta(seconds=5)
