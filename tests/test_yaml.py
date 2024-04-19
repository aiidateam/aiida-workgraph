import aiida
from aiida_workgraph import WorkGraph
import os
import pytest

aiida.load_profile()

cwd = os.path.dirname(os.path.abspath(__file__))


def test_calcfunction():
    wg = WorkGraph.from_yaml(os.path.join(cwd, "datas/test_calcfunction.yaml"))
    wg.submit(wait=True)
    assert wg.nodes["sumdiff2"].node.outputs.sum == 9


# skip this test for now
@pytest.mark.skip(reason="need to fix the identifier for a node from build_node")
def test_calcjob():
    wg = WorkGraph.from_yaml(os.path.join(cwd, "datas/test_calcjob.yaml"))
    wg.submit(wait=True)
    assert wg.nodes["add2"].node.outputs.sum == 9
