from aiida_workgraph import WorkGraph
import os
import pytest

cwd = os.path.dirname(os.path.abspath(__file__))


def test_calcfunction():
    wg = WorkGraph.from_yaml(os.path.join(cwd, "datas/test_calcfunction.yaml"))
    assert wg.tasks.float1.inputs.value.socket_value == 3.0
    assert wg.tasks.sumdiff1.inputs.x.socket_value == 2.0
    assert wg.tasks.sumdiff2.inputs.x.socket_value == 4.0
    wg.run()
    assert wg.tasks.sumdiff2.node.outputs.sum == 9


# skip this test for now
@pytest.mark.skip(reason="need to fix the identifier for a node from build_task")
def test_calcjob():
    wg = WorkGraph.from_yaml(os.path.join(cwd, "datas/test_calcjob.yaml"))
    wg.submit(wait=True)
    assert wg.tasks.add2.node.outputs.sum == 9
