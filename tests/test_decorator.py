import pytest
from aiida_workgraph import WorkGraph
from typing import Callable


def test_args() -> None:
    from aiida_workgraph import task

    @task.calcfunction()
    def test(a, b=1, **c):
        print(a, b, c)

    metadata_kwargs = set(
        [
            f"metadata.{key}"
            for key in test.process_class.spec().inputs.ports["metadata"].ports.keys()
        ]
    )
    kwargs = set(test.process_class.spec().inputs.ports.keys()).union(metadata_kwargs)
    kwargs.remove("a")
    #
    n = test.task()
    assert n.args == ["a"]
    assert set(n.kwargs) == set(kwargs)
    assert n.var_args is None
    assert n.var_kwargs == "c"
    assert n.outputs.keys() == ["result", "_outputs", "_wait"]


def test_inputs_outputs_workchain() -> None:
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    wg = WorkGraph()
    task = wg.add_task(MultiplyAddWorkChain)
    assert "metadata" in task.inputs.keys()
    assert "metadata.call_link_label" in task.inputs.keys()
    assert "result" in task.outputs.keys()


@pytest.mark.usefixtures("started_daemon_client")
def test_decorator_calcfunction(decorated_add: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_calcfunction")
    wg.add_task(decorated_add, "add1", x=2, y=3)
    wg.submit(wait=True, timeout=100)
    assert wg.tasks["add1"].outputs["result"].value == 5


def test_decorator_workfunction(decorated_add_multiply: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name="test_decorator_workfunction")
    wg.add_task(decorated_add_multiply, "add_multiply1", x=2, y=3, z=4)
    wg.submit(wait=True, timeout=100)
    assert wg.tasks["add_multiply1"].outputs["result"].value == 20


@pytest.mark.usefixtures("started_daemon_client")
def test_decorator_graph_builder(decorated_add_multiply_group: Callable) -> None:
    """Test graph build."""
    wg = WorkGraph("test_graph_builder")
    add1 = wg.add_task("AiiDAAdd", "add1", x=2, y=3)
    add_multiply1 = wg.add_task(decorated_add_multiply_group, "add_multiply1", y=3, z=4)
    sum_diff1 = wg.add_task("AiiDASumDiff", "sum_diff1")
    wg.add_link(add1.outputs[0], add_multiply1.inputs["x"])
    wg.add_link(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    wg.submit(wait=True)
    assert wg.tasks["add_multiply1"].process.outputs.result.value == 32
    assert wg.tasks["add_multiply1"].outputs["result"].value == 32
    assert wg.tasks["sum_diff1"].outputs["sum"].value == 32
