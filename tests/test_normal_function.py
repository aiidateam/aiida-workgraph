import pytest
from aiida_workgraph import WorkGraph
from typing import Callable


def test_normal_function_run(
    decorated_normal_add: Callable, decorated_add: Callable
) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_run")
    add1 = wg.tasks.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.tasks.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs["result"], add2.inputs["y"])
    wg.run()
    assert wg.tasks["add2"].node.outputs.result == 11


@pytest.mark.usefixtures("started_daemon_client")
def test_normal_function_submit(
    decorated_normal_add: Callable, decorated_add: Callable
) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name="test_normal_function_submit")
    add1 = wg.tasks.new(decorated_normal_add, "add1", x=2, y=3)
    add2 = wg.tasks.new(decorated_add, "add2", x=6)
    wg.links.new(add1.outputs["result"], add2.inputs["y"])
    wg.submit(wait=True)
    assert wg.tasks["add2"].node.outputs.result == 11
