import aiida
from utils import wait_nt

aiida.load_profile()


def test_decorator_calcfunction(decorated_add):
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    nt = WorkTree(name="test_decorator_calcfunction")
    nt.nodes.new(decorated_add, "add1", x=2, y=3)
    nt.submit(wait=True, timeout=100)
    # print("results: ", results[])
    assert nt.nodes["add1"].node.outputs.result == 5


def test_decorator_workfunction(decorated_add_multiply):
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    nt = WorkTree(name="test_decorator_workfunction")
    nt.nodes.new(decorated_add_multiply, "add_multiply1", x=2, y=3, z=4)
    nt.submit()
    wait_nt(nt, timeout=100)
    assert nt.nodes["add_multiply1"].node.outputs.result == 20


def test_decorator_node_group(decorated_add_multiply_group):
    from aiida_worktree import WorkTree

    nt = WorkTree("test_node_group")
    add1 = nt.nodes.new("AiiDAAdd", "add1", x=2, y=3)
    add_multiply1 = nt.nodes.new(
        decorated_add_multiply_group, "add_multiply1", y=3, z=4
    )
    sum_diff1 = nt.nodes.new("AiiDASumDiff", "sum_diff1")
    nt.links.new(add1.outputs[0], add_multiply1.inputs["x"])
    nt.links.new(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    nt.submit(wait=True)
    assert nt.nodes["add_multiply1"].node.outputs.group_outputs.result == 32
    assert nt.nodes["sum_diff1"].node.outputs.sum == 32
