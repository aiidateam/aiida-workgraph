import pytest
from aiida_workgraph import WorkGraph


@pytest.mark.usefixtures("started_daemon_client")
def test_build_task_from_workgraph(wg_calcfunction, decorated_add):

    wg = WorkGraph("build_task_from_workgraph")
    add1_task = wg.add_task(decorated_add, name="add1", x=1, y=3)
    wg_task = wg.add_task(wg_calcfunction, name="wg_calcfunction")
    wg.add_task(decorated_add, name="add2", y=3)
    wg.add_link(add1_task.outputs["result"], wg_task.inputs["sumdiff1.x"])
    wg.add_link(wg_task.outputs["sumdiff2.sum"], wg.tasks["add2"].inputs["x"])
    assert len(wg_task.inputs) == 7
    assert len(wg_task.outputs) == 8
    wg.submit(wait=True)
    assert wg.tasks["add2"].outputs["result"].value.value == 14
