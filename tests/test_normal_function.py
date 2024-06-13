import aiida
from aiida_workgraph import WorkGraph
from typing import Callable

aiida.load_profile()


def test_normal_function_run(
    decorated_normal_add: Callable, decorated_add: Callable
) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_run")
    add1 = wg.tasks.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.tasks.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs["result"], add2.inputs["y"])
    wg.run()
    assert wg.nodes["add2"].node.outputs.result == 11


def test_normal_function_submit(
    decorated_normal_add: Callable, decorated_add: Callable
) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_submit")
    add1 = wg.tasks.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.tasks.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs["result"], add2.inputs["y"])
    wg.submit(wait=True)
    assert wg.nodes["add2"].node.outputs.result == 11
