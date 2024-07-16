import pytest
from typing import Callable


@pytest.mark.usefixtures("started_daemon_client")
def test_node_wait(decorated_add: Callable) -> None:
    """Run simple calcfunction."""
    from aiida_workgraph import WorkGraph

    wg = WorkGraph(name="test_node_wait")
    add1 = wg.add_task(decorated_add, "add1", x=1, y=1)
    add1.set_context({"result": "sum1"})
    add2 = wg.add_task(decorated_add, "add2", x=2, y=2)
    add2.set_context({"result": "sum2"})
    add3 = wg.add_task(decorated_add, "add3", x="{{sum1}}", y="{{sum2}}")
    add3.wait = ["add1", add2]
    wg.submit(wait=True)
    add3.ctime < add1.ctime
    add3.ctime < add2.ctime
    assert add3.outputs[0].value == 6
