import pytest
from aiida_workgraph import WorkGraph, task
from typing import Callable
from aiida import orm
from node_graph import spec


def test_normal_task(decorated_add) -> None:
    """Test a normal task."""

    @task(outputs=["sum", "diff"])
    def sum_diff(x, y):
        return x + y, x - y

    wg = WorkGraph("test_normal_task")
    task1 = wg.add_task(sum_diff, name="sum_diff", x=2, y=3)
    task2 = wg.add_task(
        decorated_add, name="add", x=task1.outputs.sum, y=task1.outputs["diff"]
    )
    wg.run()
    print("node: ", task2.outputs.result)
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
def test_task_wait(decorated_add: Callable, capsys) -> None:
    """Run a WorkGraph with a task that waits on other tasks."""

    wg = WorkGraph(name="test_task_wait")
    add1 = wg.add_task(decorated_add, "add1", x=1, y=1)
    add2 = wg.add_task(decorated_add, "add2", x=2, y=2)
    add2.waiting_on.add(add1)
    wg.run()
    captured = capsys.readouterr()
    report = captured.out
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
        data["workgraph_data"]["tasks"]["add1"]["inputs"]["sockets"]["y"]["property"][
            "value"
        ]
        == 2
    )
    assert (
        data["workgraph_data"]["tasks"]["add1"]["inputs"]["sockets"]["metadata"][
            "sockets"
        ]["store_provenance"]["property"]["value"]
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
        outputs=spec.namespace(
            add_multiply=spec.namespace(add=any, multiply=any), minus=any
        )
    )
    def myfunc(x, y):
        return {
            "add_multiply": {"add": orm.Float(x + y), "multiply": orm.Float(x * y)},
            "minus": orm.Float(x - y),
        }

    wg = WorkGraph("test_namespace_outputs")
    wg.add_task(myfunc, name="myfunc", x=1.0, y=2.0)
    print(wg.tasks.myfunc.outputs)
    wg.run()
    assert wg.tasks.myfunc.outputs.minus.value == -1
    assert wg.tasks.myfunc.outputs.add_multiply.add.value == 3
    assert wg.tasks.myfunc.outputs.add_multiply.multiply.value == 2


def test_task_from_builder_add(add_code) -> None:
    """Test adding a task from a ``ProcessBuilder`` for `ArithmeticAdd`."""
    from aiida_workgraph.sockets.builtins import SocketAny
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    # Test with the ArithmeticAddCalculation

    add_builder = ArithmeticAddCalculation.get_builder()

    add_builder.code = add_code
    add_builder.x = orm.Int(2)
    add_builder.y = orm.Int(3)

    add_wg = WorkGraph("add-builder")
    add_task_name = "add_from_builder"
    add_task = add_wg.add_task(add_builder, name=add_task_name)

    assert add_task.name == add_task_name
    assert add_task.identifier == "ArithmeticAddCalculation"
    assert isinstance(add_task.inputs.x, SocketAny)
    assert add_task.inputs.x.value == 2
    assert add_wg.tasks[add_task_name].inputs["y"].value == 3
    assert add_wg.tasks[add_task_name].inputs["code"].value == add_code


def test_task_from_builder_multiply_add(add_code, decorated_add) -> None:
    """Test adding a task from a ``ProcessBuilder`` for ``MultiplyAdd``."""
    from aiida_workgraph.sockets.builtins import SocketAny
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    multiply_add_builder = MultiplyAddWorkChain.get_builder()

    multiply_add_builder.code = orm.load_code("add@localhost")
    multiply_add_builder.x = orm.Int(2)
    multiply_add_builder.y = orm.Int(3)
    multiply_add_builder.z = orm.Int(4)

    wg = WorkGraph("multiply-add-builder")

    multiply_add_task_name = "multiply_add_builder"
    add_task_name = "decorated_add"

    multiply_add_task = wg.add_task(multiply_add_builder, name=multiply_add_task_name)

    assert multiply_add_task.name == multiply_add_task_name
    assert multiply_add_task.identifier == "MultiplyAddWorkChain"
    assert multiply_add_task.inputs.x.value == 2
    assert wg.tasks[multiply_add_task_name].inputs.y.value == 3
    assert wg.tasks[multiply_add_task_name].inputs["z"].value == 4
    assert wg.tasks[multiply_add_task_name].inputs["code"].value == add_code

    # Check if task coming from ProcessBuilder behaves as other Tasks
    add_task = wg.add_task(
        decorated_add, name=add_task_name, x=multiply_add_task.outputs.result, y=11
    )

    assert isinstance(add_task.inputs.x, SocketAny)
    assert add_task.inputs.x.value is None

    assert len(wg.tasks) == 2
    assert len(wg.links) == 1
    assert wg.links_to_dict() == [
        {
            "from_node": multiply_add_task_name,
            "from_socket": "result",
            "to_node": add_task_name,
            "to_socket": "x",
        }
    ]

    wg.run()

    # Top-level test for the expected result
    assert wg.tasks.decorated_add.outputs.result.value.value == 21


def test_task_children():
    wg = WorkGraph()
    zone1 = wg.add_task("workgraph.zone", "zone1")
    zone2 = wg.add_task("workgraph.zone", "zone2")
    zone3 = wg.add_task("workgraph.zone", "zone3")
    zone2.children.add(zone1)
    with pytest.raises(ValueError, match="Task is already a child of the task: "):
        zone3.children.add(zone1)
