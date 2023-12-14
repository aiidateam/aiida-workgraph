import aiida

aiida.load_profile()


def test_args():
    from aiida_worktree import node

    @node.calcfunction()
    def test(a, b=1, **c):
        print(a, b, c)

    #
    n = test.node()
    assert n.args == []
    assert n.kwargs == ["a", "b"]
    assert n.var_args is None
    assert n.var_kwargs == "c"


def test_decorator_calcfunction(decorated_add):
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    wt = WorkTree(name="test_decorator_calcfunction")
    wt.nodes.new(decorated_add, "add1", x=2, y=3)
    wt.submit(wait=True, timeout=100)
    assert wt.nodes["add1"].node.outputs.result == 5


def test_decorator_workfunction(decorated_add_multiply):
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    wt = WorkTree(name="test_decorator_workfunction")
    wt.nodes.new(decorated_add_multiply, "add_multiply1", x=2, y=3, z=4)
    wt.submit(wait=True, timeout=100)
    assert wt.nodes["add_multiply1"].node.outputs.result == 20


def test_decorator_node_group(decorated_add_multiply_group):
    from aiida_worktree import WorkTree

    wt = WorkTree("test_node_group")
    add1 = wt.nodes.new("AiiDAAdd", "add1", x=2, y=3, t=10)
    add_multiply1 = wt.nodes.new(
        decorated_add_multiply_group, "add_multiply1", y=3, z=4
    )
    sum_diff1 = wt.nodes.new("AiiDASumDiff", "sum_diff1")
    wt.links.new(add1.outputs[0], add_multiply1.inputs["x"])
    wt.links.new(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    wt.submit(wait=True)
    assert wt.nodes["add_multiply1"].node.outputs.group_outputs.result == 32
    assert wt.nodes["sum_diff1"].node.outputs.sum == 32
