from aiida_workgraph import WorkGraph
from aiida.cmdline.utils.common import get_workchain_report


def test_zone_task(decorated_add, capsys):
    """Test the zone task."""

    wg = WorkGraph("test_zone")
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    zone1 = wg.add_task("workgraph.zone", name="zone1")
    zone1.add_task(decorated_add, name="add2", x=1, y=1)
    zone1.add_task(decorated_add, name="add3", x=1, y=add1.outputs.result)
    wg.add_task(decorated_add, name="add4", x=1, y=wg.tasks.add2.outputs.result)
    wg.add_task(decorated_add, name="add5", x=1, y=wg.tasks.add3.outputs.result)
    wg.run()
    captured = capsys.readouterr()
    report = captured.out
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add2,add3" in report
    assert "tasks ready to run: add4,add5" in report
    # load the WorkGraph should add the cihld tasks
    wg = WorkGraph.load(wg.process.pk)
    assert len(wg.tasks.zone1.children) == 2
