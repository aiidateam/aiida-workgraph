import pytest
from typing import Callable
from aiida_workgraph import WorkGraph
from aiida.cmdline.utils.common import get_workchain_report
from aiida_workgraph.engine.scheduler.client import get_scheduler
from aiida import orm


@pytest.mark.skip("Skip for now")
@pytest.mark.usefixtures("started_daemon_client")
def test_scheduler(decorated_add: Callable, started_scheduler_client) -> None:
    """Test graph build."""
    wg = WorkGraph("test_scheduler")
    add1 = wg.add_task(decorated_add, x=2, y=3)
    add2 = wg.add_task(decorated_add, "add2", x=3, y=add1.outputs["result"])
    # use run to check if graph builder workgraph can be submit inside the engine
    wg.submit(to_scheduler=True, wait=True)
    pk = get_scheduler()
    report = get_workchain_report(orm.load(pk), "REPORT")
    print("report: ", report)
    assert add2.outputs["result"].value == 8
