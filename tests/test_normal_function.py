import aiida
from aiida_workgraph import WorkGraph

aiida.load_profile()


def test_normal_function_run(decorated_normal_add, decorated_add):
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_run")
    add1 = wg.nodes.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.nodes.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs[0], add2.inputs[1])
    wg.run()
    assert wg.nodes["add2"].node.outputs.result == 11


def test_normal_function_submit(decorated_normal_add, decorated_add):
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_submit")
    add1 = wg.nodes.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.nodes.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs[0], add2.inputs[1])
    wg.submit(wait=True)
    assert wg.nodes["add2"].node.outputs.result == 11
