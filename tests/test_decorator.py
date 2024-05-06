import aiida
from aiida_workgraph import WorkGraph
from typing import Callable

aiida.load_profile()


def test_args() -> None:
    from aiida_workgraph import node

    @node.calcfunction()
    def test(a, b=1, **c):
        print(a, b, c)

    #
    n = test.node()
    assert n.args == ["a"]
    assert n.kwargs == ["b"]
    assert n.var_args is None
    assert n.var_kwargs == "c"


def test_decorator_calcfunction(decorated_add: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_calcfunction")
    wg.nodes.new(decorated_add, "add1", x=2, y=3)
    wg.submit(wait=True, timeout=100)
    assert wg.nodes["add1"].outputs["result"].value == 5


def test_decorator_workfunction(decorated_add_multiply: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_workfunction")
    wg.nodes.new(decorated_add_multiply, "add_multiply1", x=2, y=3, z=4)
    wg.submit(wait=True, timeout=100)
    assert wg.nodes["add_multiply1"].outputs["result"].value == 20


def test_decorator_graph_builder(decorated_add_multiply_group: Callable) -> None:
    """Test graph build."""
    wg = WorkGraph("test_graph_builder")
    add1 = wg.nodes.new("AiiDAAdd", "add1", x=2, y=3, t=10)
    add_multiply1 = wg.nodes.new(
        decorated_add_multiply_group, "add_multiply1", y=3, z=4
    )
    sum_diff1 = wg.nodes.new("AiiDASumDiff", "sum_diff1")
    wg.links.new(add1.outputs[0], add_multiply1.inputs["x"])
    wg.links.new(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    wg.submit(wait=True)
    assert wg.nodes["add_multiply1"].outputs["result"].value == 32
    assert wg.nodes["sum_diff1"].outputs["sum"].value == 32
