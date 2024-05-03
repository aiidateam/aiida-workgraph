from aiida_workgraph import WorkGraph, build_node
from aiida import load_profile, orm
import time
import pytest

load_profile()


def test_to_dict(wg_calcjob):
    """Export NodeGraph to dict."""
    wg = wg_calcjob
    ntdata = wg.to_dict()
    assert len(ntdata["nodes"]) == len(wg.nodes)
    assert len(ntdata["links"]) == len(wg.links)


def test_from_dict(wg_calcjob):
    """Export NodeGraph to dict."""
    wg = wg_calcjob
    ntdata = wg.to_dict()
    wg1 = WorkGraph.from_dict(ntdata)
    assert len(wg.nodes) == len(wg1.nodes)
    assert len(wg.links) == len(wg1.links)


def test_new_node(wg_calcjob, arithmetic_add):
    """Add new node."""
    wg = wg_calcjob
    n = len(wg.nodes)
    wg.nodes.new(arithmetic_add)
    assert len(wg.nodes) == n + 1


def test_save_load(wg_calcjob):
    """Save the workgraph"""
    wg = wg_calcjob
    wg.save()
    assert wg.process.process_state.value.upper() == "CREATED"
    wg2 = WorkGraph.load(wg.process.pk)
    assert len(wg.nodes) == len(wg2.nodes)


# skip this test
@pytest.mark.skip(reason="PAUSED state is wrong for the moment.")
def test_pause(wg_engine):
    wg = wg_engine
    wg.name = "test_pause"
    wg.submit()
    time.sleep(5)
    wg.pause()
    wg.update()
    assert wg.process.process_state.value.upper() == "PAUSED"


def test_reset_message(wg_calcjob):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    wg = wg_calcjob
    wg.submit()
    wg = WorkGraph.load(wg.process.pk)
    wg.nodes["add2"].set({"y": orm.Int(10).store()})
    wg.save()
    msgs = wg.process.base.extras.get("workgraph_queue", [])
    assert len(msgs) == 1


def test_restart(wg_calcjob):
    """Restart from a finished workgraph.
    Load the workgraph, modify the node, and restart the workgraph.
    Only the modified node and its child nodes will be rerun."""
    wg = wg_calcjob
    wg.name = "test_restart_0"
    wg.submit(wait=True)
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.name = "test_restart_1"
    wg1.nodes["add2"].set({"y": orm.Int(10).store()})
    wg1.submit(wait=True, restart=True)
    wg1.update()
    assert wg1.nodes["add3"].node.outputs.sum == 13
    assert wg1.nodes["add1"].node.pk == wg.nodes["add1"].pk
    assert wg1.nodes["add2"].node.pk != wg.nodes["add2"].pk


def test_append_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph("test_graph_build")
    add1 = wg.nodes.new("AiiDAAdd", "add1", x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    # extend workgraph
    wg.extend(add_multiply_wg, prefix="group_")
    wg.links.new(add1.outputs[0], wg.nodes["group_add1"].inputs["x"])
    wg.submit(wait=True)
    assert wg.nodes["group_multiply1"].node.outputs.result == 45


def test_node_from_workgraph(decorated_add_multiply_group):
    wg = WorkGraph("test_node_from_workgraph")
    add1 = wg.nodes.new("AiiDAAdd", "add1", x=2, y=3)
    add2 = wg.nodes.new("AiiDAAdd", "add2", y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    AddMultiplyNode = build_node(add_multiply_wg)
    assert "add1.x" in AddMultiplyNode().inputs.keys()
    # add the workgraph as a node
    add_multiply1 = wg.nodes.new(AddMultiplyNode, "add_multiply1")
    wg.links.new(add1.outputs[0], add_multiply1.inputs["add1.x"])
    wg.links.new(add_multiply1.outputs["multiply1.result"], add2.inputs["x"])
    # wg.submit(wait=True)
    wg.run()
    assert wg.nodes["add2"].node.outputs.sum == 48
