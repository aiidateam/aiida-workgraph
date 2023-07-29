import aiida

aiida.load_profile()


def test_worktree_ctx(decorated_add):
    """Set/get data to/from context."""
    from aiida_worktree import WorkTree
    from aiida.orm import Float

    nt = WorkTree(name="test_worktree_ctx")
    nt.ctx = {"x": Float(2)}
    add1 = nt.nodes.new(decorated_add, "add1", x="{{x}}", y=Float(3).store())
    nt.submit(wait=True)
    assert add1.node.outputs.result.value == 5


def test_node_to_ctx(decorated_add):
    """Set/get data to/from context."""
    from aiida_worktree import WorkTree
    from aiida.orm import Float

    nt = WorkTree(name="test_node_to_ctx")
    add1 = nt.nodes.new(decorated_add, "add1", x=Float(2).store(), y=Float(3).store())
    add1.to_ctx = [["result", "sum"]]
    add2 = nt.nodes.new(decorated_add, "add2", y="{{ sum }}")
    nt.links.new(add1.outputs[0], add2.inputs["x"])
    nt.submit(wait=True)
    assert add2.node.outputs.result.value == 10
