import aiida
from aiida_worktree import WorkTree

aiida.load_profile()


def test_normal_function_run(decorated_normal_add, decorated_add):
    """Run simple calcfunction."""
    nt = WorkTree(name="test_normal_function_run")
    add1 = nt.nodes.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = nt.nodes.new(decorated_add, "add2", x=6)
    nt.links.new(add1.outputs[0], add2.inputs[1])
    nt.run()
    assert nt.nodes["add2"].node.outputs.result == 11


def test_normal_function_submit(decorated_normal_add, decorated_add):
    """Run simple calcfunction."""
    nt = WorkTree(name="test_normal_function_submit")
    add1 = nt.nodes.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = nt.nodes.new(decorated_add, "add2", x=6)
    nt.links.new(add1.outputs[0], add2.inputs[1])
    nt.submit(wait=True)
    assert nt.nodes["add2"].node.outputs.result == 11
