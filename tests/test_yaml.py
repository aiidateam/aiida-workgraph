import aiida
from aiida_worktree import WorkTree
import os

aiida.load_profile()

cwd = os.path.dirname(os.path.abspath(__file__))


def test_calcfunction():
    nt = WorkTree.from_yaml(os.path.join(cwd, "datas/test_calcfunction.yaml"))
    nt.submit(wait=True)
    assert nt.nodes["sumdiff2"].node.outputs.sum == 9


def test_calcjob():
    nt = WorkTree.from_yaml(os.path.join(cwd, "datas/test_calcjob.yaml"))
    nt.submit(wait=True)
    assert nt.nodes["add2"].node.outputs.sum == 9
