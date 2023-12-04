from aiida_worktree import WorkTree
from aiida import load_profile

load_profile()


def test_to_dict(wt_calcjob):
    """Export NodeGraph to dict."""
    wt = wt_calcjob
    ntdata = wt.to_dict()
    assert len(ntdata["nodes"]) == len(wt.nodes)
    assert len(ntdata["links"]) == len(wt.links)


def test_from_dict(wt_calcjob):
    """Export NodeGraph to dict."""
    wt = wt_calcjob
    ntdata = wt.to_dict()
    wt1 = WorkTree.from_dict(ntdata)
    assert len(wt.nodes) == len(wt1.nodes)
    assert len(wt.links) == len(wt1.links)


def test_new_node(wt_calcjob, arithmetic_add):
    """Add new node."""
    wt = wt_calcjob
    n = len(wt.nodes)
    wt.nodes.new(arithmetic_add)
    assert len(wt.nodes) == n + 1


def test_save_load(wt_calcjob):
    """Save the worktree"""
    wt = wt_calcjob
    wt.save()
    assert wt.process.process_state.value.upper() == "CREATED"
    wt2 = WorkTree.load(wt.process.pk)
    assert len(wt.nodes) == len(wt2.nodes)
