from aiida_worktree import WorkTree, build_node
from aiida import load_profile, orm
import time
import pytest

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


# skip this test
@pytest.mark.skip(reason="PAUSED state is wrong for the moment.")
def test_pause(wt_engine):
    wt = wt_engine
    wt.name = "test_pause"
    wt.submit()
    time.sleep(5)
    wt.pause()
    wt.update()
    assert wt.process.process_state.value.upper() == "PAUSED"


def test_reset_message(wt_calcjob):
    """Modify a node and save the worktree.
    This will add a message to the worktree_queue extra field."""
    wt = wt_calcjob
    wt.submit()
    wt = WorkTree.load(wt.process.pk)
    wt.nodes["add2"].set({"y": orm.Int(10).store()})
    wt.save()
    msgs = wt.process.base.extras.get("worktree_queue", [])
    assert len(msgs) == 1


def test_restart(wt_calcjob):
    """Restart from a finished worktree.
    Load the worktree, modify the node, and restart the worktree.
    Only the modified node and its child nodes will be rerun."""
    wt = wt_calcjob
    wt.name = "test_restart_0"
    wt.submit(wait=True)
    wt1 = WorkTree.load(wt.process.pk)
    wt1.name = "test_restart_1"
    wt1.nodes["add2"].set({"y": orm.Int(10).store()})
    wt1.submit(wait=True, restart=True)
    wt1.update()
    assert wt1.nodes["add3"].node.outputs.sum == 13
    assert wt1.nodes["add1"].node.pk == wt.nodes["add1"].pk
    assert wt1.nodes["add2"].node.pk != wt.nodes["add2"].pk


def test_append_worktree(decorated_add_multiply_group):
    from aiida_worktree import WorkTree

    wt = WorkTree("test_node_group")
    add1 = wt.nodes.new("AiiDAAdd", "add1", x=2, y=3)
    add_multiply_wt = decorated_add_multiply_group(x=0, y=4, z=5)
    # append worktree
    wt.append(add_multiply_wt, prefix="group_")
    wt.links.new(add1.outputs[0], wt.nodes["group_add1"].inputs["x"])
    wt.submit(wait=True)
    assert wt.nodes["group_multiply1"].node.outputs.result == 45


def test_node_from_worktree(decorated_add_multiply_group):
    wt = WorkTree("test_node_from_worktree")
    add1 = wt.nodes.new("AiiDAAdd", "add1", x=2, y=3)
    add2 = wt.nodes.new("AiiDAAdd", "add2", y=3)
    add_multiply_wt = decorated_add_multiply_group(x=0, y=4, z=5)
    AddMultiplyNode = build_node(add_multiply_wt)
    assert "add1.x" in AddMultiplyNode().inputs.keys()
    # add the worktree as a node
    add_multiply1 = wt.nodes.new(AddMultiplyNode, "add_multiply1")
    wt.links.new(add1.outputs[0], add_multiply1.inputs["add1.x"])
    wt.links.new(add_multiply1.outputs["multiply1.result"], add2.inputs["x"])
    # wt.submit(wait=True)
    wt.run()
    assert wt.nodes["add2"].node.outputs.sum == 48
