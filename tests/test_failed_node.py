from typing import Callable
import aiida

aiida.load_profile()


def test_failed_node(decorated_sqrt: Callable, decorated_add: Callable) -> None:
    """Submit simple calcfunction."""
    from aiida_workgraph import WorkGraph
    from aiida.orm import Float

    wg = WorkGraph(name="test_failed_node")
    wg.nodes.new(decorated_add, "add1", x=Float(1), y=Float(2))
    wg.nodes.new(decorated_sqrt, "sqrt1", x=Float(-1))
    wg.submit(wait=True)
    # print("results: ", results[])
    assert wg.process.exit_status == 302
