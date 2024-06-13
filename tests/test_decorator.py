import aiida
from aiida_workgraph import WorkGraph
from typing import Callable

aiida.load_profile()


def test_args() -> None:
    from aiida_workgraph import task

    @task.calcfunction()
    def test(a, b=1, **c):
        print(a, b, c)

    #
    n = test.node()
    assert n.args == ["a"]
    assert n.kwargs == [
        "metadata",
        "metadata.store_provenance",
        "metadata.description",
        "metadata.label",
        "metadata.call_link_label",
        "b",
    ]
    assert n.var_args is None
    assert n.var_kwargs == "c"
    assert n.outputs.keys() == ["result", "_outputs", "_wait"]


def test_inputs_outputs_workchain() -> None:
    from aiida_quantumespresso.workflows.pdos import PdosWorkChain

    wg = WorkGraph()
    pdos = wg.tasks.new(PdosWorkChain)
    assert "scf" in pdos.inputs.keys()
    assert "scf.pw" in pdos.inputs.keys()
    assert "scf.pw.metadata" in pdos.inputs.keys()
    assert "dos" in pdos.outputs.keys()
    assert "dos.remote_folder" in pdos.outputs.keys()


def test_decorator_calcfunction(decorated_add: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_calcfunction")
    wg.tasks.new(decorated_add, "add1", x=2, y=3)
    wg.submit(wait=True, timeout=100)
    assert wg.nodes["add1"].outputs["result"].value == 5


def test_decorator_workfunction(decorated_add_multiply: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_workfunction")
    wg.tasks.new(decorated_add_multiply, "add_multiply1", x=2, y=3, z=4)
    wg.submit(wait=True, timeout=100)
    assert wg.nodes["add_multiply1"].outputs["result"].value == 20


def test_decorator_graph_builder(decorated_add_multiply_group: Callable) -> None:
    """Test graph build."""
    wg = WorkGraph("test_graph_builder")
    add1 = wg.tasks.new("AiiDAAdd", "add1", x=2, y=3, t=10)
    add_multiply1 = wg.tasks.new(
        decorated_add_multiply_group, "add_multiply1", y=3, z=4
    )
    sum_diff1 = wg.tasks.new("AiiDASumDiff", "sum_diff1")
    wg.links.new(add1.outputs[0], add_multiply1.inputs["x"])
    wg.links.new(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    wg.submit(wait=True)
    assert wg.nodes["add_multiply1"].outputs["result"].value == 32
    assert wg.nodes["sum_diff1"].outputs["sum"].value == 32
