import pytest
from aiida_workgraph import WorkGraph, task
from typing import Callable


def test_custom_outputs():
    """Test custom outputs."""

    @task(outputs=["sum", {"name": "product", "identifier": "workgraph.any"}])
    def add_multiply(x, y):
        return {"sum": x + y, "product": x * y}

    n = add_multiply.task()
    assert "sum" in n.outputs.keys()
    assert "product" in n.outputs.keys()


@pytest.fixture(params=["decorator_factory", "decorator"])
def task_calcfunction(request):
    if request.param == "decorator_factory":

        @task.calcfunction()
        def test(a, b=1, **c):
            print(a, b, c)

    elif request.param == "decorator":

        @task.calcfunction
        def test(a, b=1, **c):
            print(a, b, c)

    else:
        raise ValueError(f"{request.param} not supported.")
    return test


def test_decorators_calcfunction_args(task_calcfunction) -> None:
    metadata_kwargs = set(
        [
            f"metadata.{key}"
            for key in task_calcfunction.process_class.spec()
            .inputs.ports["metadata"]
            .ports.keys()
        ]
    )
    kwargs = set(task_calcfunction.process_class.spec().inputs.ports.keys()).union(
        metadata_kwargs
    )
    kwargs.remove("a")
    #
    n = task_calcfunction.task()
    assert n.args == ["a"]
    assert set(n.kwargs) == set(kwargs)
    assert n.var_args is None
    assert n.var_kwargs == "c"
    assert n.outputs.keys() == ["result", "_outputs", "_wait"]


@pytest.fixture(params=["decorator_factory", "decorator"])
def task_function(request):
    if request.param == "decorator_factory":

        @task()
        def test(a, b=1, **c):
            print(a, b, c)

    elif request.param == "decorator":

        @task
        def test(a, b=1, **c):
            print(a, b, c)

    else:
        raise ValueError(f"{request.param} not supported.")
    return test


def test_decorators_task_args(task_function):

    tdata = task_function.tdata
    assert tdata["args"] == ["a"]
    assert tdata["kwargs"] == ["b"]
    assert tdata["var_args"] is None
    assert tdata["var_kwargs"] == "c"
    assert set(tdata["outputs"].keys()) == set(["result", "_outputs", "_wait"])


@pytest.fixture(params=["decorator_factory", "decorator"])
def task_workfunction(request):
    if request.param == "decorator_factory":

        @task.workfunction()
        def test(a, b=1, **c):
            print(a, b, c)

    elif request.param == "decorator":

        @task.workfunction
        def test(a, b=1, **c):
            print(a, b, c)

    else:
        raise ValueError(f"{request.param} not supported.")
    return test


def test_decorators_workfunction_args(task_workfunction) -> None:
    metadata_kwargs = set(
        [
            f"metadata.{key}"
            for key in task_workfunction.process_class.spec()
            .inputs.ports["metadata"]
            .ports.keys()
        ]
    )
    kwargs = set(task_workfunction.process_class.spec().inputs.ports.keys()).union(
        metadata_kwargs
    )
    kwargs.remove("a")
    #
    n = task_workfunction.task()
    assert n.args == ["a"]
    assert set(n.kwargs) == set(kwargs)
    assert n.var_args is None
    assert n.var_kwargs == "c"
    assert n.outputs.keys() == ["result", "_outputs", "_wait"]


def test_decorators_parameters() -> None:
    """Test passing parameters to decorators."""

    @task.calcfunction(
        inputs=[{"name": "c", "link_limit": 1000}],
        outputs=[{"name": "sum"}, {"name": "product"}],
    )
    def test(a, b=1, **c):
        return {"sum": a + b, "product": a * b}

    test1 = test.task()
    assert test1.inputs["c"].link_limit == 1000
    assert "sum" in test1.outputs.keys()
    assert "product" in test1.outputs.keys()


@pytest.fixture(params=["decorator_factory", "decorator"])
def task_graph_builder(request):
    if request.param == "decorator_factory":

        @task.graph_builder()
        def add_multiply_group(a, b=1, **c):
            wg = WorkGraph("add_multiply_group")
            print(a, b, c)
            return wg

    elif request.param == "decorator":

        @task.graph_builder
        def add_multiply_group(a, b=1, **c):
            wg = WorkGraph("add_multiply_group")
            print(a, b, c)
            return wg

    else:
        raise ValueError(f"{request.param} not supported.")

    return add_multiply_group


def test_decorators_graph_builder_args(task_graph_builder) -> None:
    assert task_graph_builder.identifier == "add_multiply_group"
    n = task_graph_builder.task()
    assert n.args == ["a"]
    assert n.kwargs == ["b"]
    assert n.var_args is None
    assert n.var_kwargs == "c"
    assert set(n.outputs.keys()) == set(["_outputs", "_wait"])


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
    add1 = wg.add_task("workgraph.test_add", "add1", x=2, y=3)
    add_multiply1 = wg.add_task(decorated_add_multiply_group, "add_multiply1", y=3, z=4)
    sum_diff1 = wg.add_task("workgraph.test_sum_diff", "sum_diff1")
    wg.add_link(add1.outputs[0], add_multiply1.inputs["x"])
    wg.add_link(add_multiply1.outputs["result"], sum_diff1.inputs["x"])
    # use run to check if graph builder workgraph can be submit inside the engine
    wg.run()
    assert wg.tasks["add_multiply1"].process.outputs.result.value == 32
    assert wg.tasks["add_multiply1"].outputs["result"].value == 32
    assert wg.tasks["sum_diff1"].outputs["sum"].value == 32
