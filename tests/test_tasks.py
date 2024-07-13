import aiida
from aiida_workgraph import WorkGraph


aiida.load_profile()


def test_build_task_from_workgraph(wg_calcfunction, decorated_add):

    wg = WorkGraph("build_task_from_workgraph")
    add1_task = wg.tasks.new(decorated_add, name="add1", x=1, y=3)
    wg_task = wg.tasks.new(wg_calcfunction, name="wg_calcfunction")
    wg.tasks.new(decorated_add, name="add2", y=3)
    wg.links.new(add1_task.outputs["result"], wg_task.inputs["sumdiff1.x"])
    wg.links.new(wg_task.outputs["sumdiff3.sum"], wg.tasks["add2"].inputs["x"])
    assert len(wg_task.inputs) == 15
    assert len(wg_task.outputs) == 13
    wg.submit(wait=True)
    assert wg.tasks["add2"].outputs["result"].value.value == 20
