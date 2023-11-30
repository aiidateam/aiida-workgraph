import aiida
from aiida_worktree import WorkTree
import os
import pytest

aiida.load_profile()

cwd = os.path.dirname(os.path.abspath(__file__))


def test_calcfunction():
    wt = WorkTree.from_yaml(os.path.join(cwd, "datas/test_calcfunction.yaml"))
    wt.submit(wait=True)
    assert wt.nodes["sumdiff2"].node.outputs.sum == 9


# skip this test for now
@pytest.mark.skip(reason="need to fix the identifier for a node from build_node")
def test_calcjob():
    wt = WorkTree.from_yaml(os.path.join(cwd, "datas/test_calcjob.yaml"))
    wt.submit(wait=True)
    assert wt.nodes["add2"].node.outputs.sum == 9
