import pytest
from typing import Callable


@pytest.mark.usefixtures("started_daemon_client")
def test_failed_node(decorated_sqrt: Callable, decorated_add: Callable) -> None:
    """Submit simple calcfunction."""
    from aiida_workgraph import WorkGraph
    from aiida.orm import Float

    wg = WorkGraph(name="test_failed_node")
    wg.tasks.new(decorated_add, "add1", x=Float(1), y=Float(2))
    sqrt1 = wg.tasks.new(decorated_sqrt, "sqrt1", x=Float(-1))
    wg.tasks.new(decorated_sqrt, "sqrt2", x=sqrt1.outputs["result"])
    wg.submit(wait=True)
    # print("results: ", results[])
    assert wg.process.exit_status == 302
    assert (
        wg.process.exit_message
        == "WorkGraph finished, but tasks: ['sqrt1'] failed. Thus all their child tasks are skipped."
    )
