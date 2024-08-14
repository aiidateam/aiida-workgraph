from aiida_workgraph import WorkGraph
from aiida.cmdline.utils.common import get_workchain_report


def test_zone_task(decorated_add):
    """Test the zone task."""

    wg = WorkGraph("test_zone")
    wg.context = {}
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    wg.add_task(decorated_add, name="add2", x=1, y=1)
    add3 = wg.add_task(decorated_add, name="add3", x=1, y=add1.outputs["result"])
    wg.add_task(decorated_add, name="add4", x=1, y=add3.outputs["result"])
    wg.add_task(decorated_add, name="add5", x=1, y=add3.outputs["result"])
    zone1 = wg.add_task("workgraph.zone", name="Zone1")
    zone1.children = ["add2", "add3"]
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add1" in report
    assert "tasks ready to run: add2,add3" in report
    assert "tasks ready to run: add4" in report
    assert "tasks ready to run: add5" in report
