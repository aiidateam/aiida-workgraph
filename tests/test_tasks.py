import pytest
from aiida_workgraph import WorkGraph, task
from typing import Callable
from aiida.cmdline.utils.common import get_workchain_report
from aiida import orm


def test_normal_task(decorated_add) -> None:
    """Test a normal task."""

    @task(outputs=[{"name": "sum"}, {"name": "diff"}])
    def sum_diff(x, y):
        return x + y, x - y

    wg = WorkGraph("test_normal_task")
    task1 = wg.add_task(sum_diff, name="sum_diff", x=2, y=3)
    task2 = wg.add_task(
        decorated_add, name="add", x=task1.outputs.sum, y=task1.outputs["diff"]
    )
    wg.run()
    print("node: ", task2.node.outputs.result)
    wg.update()
    assert task2.outputs.result.value == 4


def test_task_collection(decorated_add: Callable) -> None:
    """Test the TaskCollection class.
    Since waiting_on and children are TaskCollection instances, we test the waiting_on."""
    wg = WorkGraph("test_task_collection")
    for i in range(5):
        wg.add_task(decorated_add, name=f"task{i}")
    task1 = wg.tasks["task1"]
    # check the graph is not None
    assert task1.waiting_on is not None
    # add a task to waiting_on
    task1.waiting_on.add("task2")
    task1.waiting_on.add(["task3", wg.tasks["task4"]])
    assert len(task1.waiting_on) == 3
    assert wg.tasks["task4"] in task1.waiting_on
    # remove a task from waiting_on
    task1.waiting_on.remove("task2")
    assert wg.tasks["task2"] not in task1.waiting_on
    # clear waiting_on
    task1.waiting_on.clear()
    assert len(task1.waiting_on) == 0


@pytest.mark.usefixtures("started_daemon_client")
def test_task_wait(decorated_add: Callable) -> None:
    """Run a WorkGraph with a task that waits on other tasks."""

    wg = WorkGraph(name="test_task_wait")
    add1 = wg.add_task(decorated_add, "add1", x=1, y=1)
    add2 = wg.add_task(decorated_add, "add2", x=2, y=2)
    add2.waiting_on.add(add1)
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "tasks ready to run: add1" in report


def test_set_non_dynamic_namespace_socket(decorated_add) -> None:
    """Test setting the namespace of a task."""
    from .utils.test_workchain import WorkChainWithNestNamespace

    wg = WorkGraph(name="test_set_namespace")
    task1 = wg.add_task(decorated_add)
    task2 = wg.add_task(WorkChainWithNestNamespace)
    task2.set(
        {
            "non_dynamic_port": {"a": task1.outputs.result, "b": orm.Int(2)},
        }
    )
    assert len(task2.inputs["non_dynamic_port.a"]._links) == 1
    assert task2.inputs["non_dynamic_port"]._value == {"b": orm.Int(2)}
    assert len(wg.links) == 1


def test_set_namespace_socket(decorated_add) -> None:
    """Test setting the namespace of a task."""
    from .utils.test_workchain import WorkChainWithNestNamespace

    wg = WorkGraph(name="test_set_namespace")
    task1 = wg.add_task(decorated_add)
    task2 = wg.add_task(WorkChainWithNestNamespace)
    task2.set(
        {
            "add": {"x": task1.outputs.result, "y": orm.Int(2)},
        }
    )
    assert len(task2.inputs["add.x"]._links) == 1
    assert task2.inputs["add"]._value == {
        "y": orm.Int(2),
    }
    assert len(wg.links) == 1


def test_set_dynamic_port_input(decorated_add) -> None:
    """Test setting dynamic port input of a task.
    Use can pass AiiDA nodes as values of the dynamic port,
    and the task will create the input for each item in the dynamic port.
    """
    from .utils.test_workchain import WorkChainWithNestNamespace

    wg = WorkGraph(name="test_set_dynamic_port_input")
    task1 = wg.add_task(decorated_add)
    task2 = wg.add_task(
        WorkChainWithNestNamespace,
        dynamic_port={
            "input1": None,
            "input2": orm.Int(2),
            "input3": task1.outputs.result,
            "nested": {"input4": orm.Int(4), "input5": task1.outputs.result},
        },
    )
    wg.add_link(task1.outputs["_wait"], task2.inputs["dynamic_port.input1"])
    # task will create input for each item in the dynamic port (nodes)
    assert "dynamic_port.input1" in task2.inputs
    assert "dynamic_port.input2" in task2.inputs
    # if the value of the item is a Socket, then it will create a link, and pop the item
    assert "dynamic_port.input3" in task2.inputs
    assert "dynamic_port.nested.input4" in task2.inputs
    assert "dynamic_port.nested.input5" in task2.inputs
    assert task2.inputs.dynamic_port._value == {
        "input2": orm.Int(2),
        "nested": {"input4": orm.Int(4)},
    }
    assert len(wg.links) == 3


def test_set_dynamic_port_output(add_code) -> None:
    """Test setting dynamic port input of a task.
    Use can pass AiiDA nodes as values of the dynamic port,
    and the task will create the input for each item in the dynamic port.
    """
    from .utils.test_workchain import WorkChainWithNestNamespace

    wg = WorkGraph(name="test_set_dynamic_port_input")
    task1 = wg.add_task(WorkChainWithNestNamespace, name="task1")
    task1.set(
        {
            "add": {"x": orm.Int(1), "y": orm.Int(2), "code": add_code},
            "multiply_add": {
                "x": orm.Int(1),
                "y": orm.Int(2),
                "z": orm.Int(3),
                "code": add_code,
            },
        }
    )
    wg.run()
    assert wg.tasks.task1.outputs.dynamic_port.result1.value == 5


def test_set_inputs(decorated_add: Callable) -> None:
    """Test setting inputs of a task."""

    wg = WorkGraph(name="test_set_inputs")
    add1 = wg.add_task(decorated_add, "add1", x=1)
    add1.set({"y": 2, "metadata.store_provenance": False})
    data = wg.prepare_inputs(metadata=None)
    assert (
        data["workgraph_data"]["tasks"]["add1"]["inputs"]["y"]["property"]["value"] == 2
    )
    assert (
        data["workgraph_data"]["tasks"]["add1"]["inputs"]["metadata"]["sockets"][
            "store_provenance"
        ]["property"]["value"]
        is False
    )


def test_set_inputs_from_builder(add_code) -> None:
    """Test setting inputs of a task from a builder function."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    wg = WorkGraph(name="test_set_inputs_from_builder")
    add1 = wg.add_task(ArithmeticAddCalculation, "add1")
    # create the builder
    builder = add_code.get_builder()
    builder.x = 1
    builder.y = 2
    add1.set_from_builder(builder)
    assert add1.inputs.x.value == 1
    assert add1.inputs["y"].value == 2
    assert add1.inputs["code"].value == add_code
    with pytest.raises(
        AttributeError,
        match=f"Executor {ArithmeticAddCalculation.__name__} does not have the get_builder_from_protocol method.",
    ):
        add1.set_from_protocol(code=add_code, protocol="fast")


def test_namespace_outputs():
    @task.calcfunction(
        outputs=[
            {"identifier": "workgraph.namespace", "name": "add_multiply"},
            {"name": "add_multiply.add"},
            {"name": "add_multiply.multiply"},
            {"name": "minus"},
        ]
    )
    def myfunc(x, y):
        return {
            "add_multiply": {"add": orm.Float(x + y), "multiply": orm.Float(x * y)},
            "minus": orm.Float(x - y),
        }

    wg = WorkGraph("test_namespace_outputs")
    wg.add_task(myfunc, name="myfunc", x=1.0, y=2.0)
    wg.run()
    assert wg.tasks.myfunc.outputs.minus.value == -1
    assert wg.tasks.myfunc.outputs.add_multiply.add.value == 3
    assert wg.tasks.myfunc.outputs.add_multiply.multiply.value == 2
