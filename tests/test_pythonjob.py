import pytest
from aiida_workgraph import WorkGraph, task, Task
from typing import Any
import numpy as np
from aiida_workgraph.socket_spec import namespace


def test_to_dict():
    """Test the to_dict method.
    PythonJobTask override the `to_dict` method to serialize the raw
    data to AiiDA data."""
    from aiida import orm

    @task.pythonjob()
    def repeat_data(data: list, dim: int = 2) -> list:
        return data * dim

    wg = WorkGraph("test_to_dict")
    # atoms will be converted to AtomsData automatically
    wg.add_task(
        repeat_data,
        data=[1, 2, 3],
        dim=2,
    )
    data = wg.tasks.repeat_data.to_dict()
    assert not isinstance(
        data["inputs"]["data"],
        orm.Data,
    )
    data = wg.tasks.repeat_data.to_dict(should_serialize=True)
    assert isinstance(
        data["inputs"]["data"],
        orm.Data,
    )


def test_imported_pythonjob(fixture_localhost, python_executable_path):
    from aiida_workgraph.executors.test import add_pythonjob
    from aiida import orm

    wg = WorkGraph("test_imported_pythonjob")
    wg.add_task(
        add_pythonjob,
        name="add1",
        x=1,
        y=2,
        computer="localhost",
        command_info={"label": python_executable_path},
    )
    wg.run()
    assert isinstance(wg.tasks.add1.outputs.result.value, orm.Data)
    assert wg.tasks.add1.outputs.result.value == 3


@pytest.mark.usefixtures("started_daemon_client")
def test_decorator(fixture_localhost, python_executable_path):
    """Test decorator."""

    @task.pythonjob(
        outputs=["sum", "diff"],
    )
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    def multiply(x: Any, y: Any) -> Any:
        return x * y

    decorted_multiply = task.pythonjob()(multiply)

    wg = WorkGraph("test_PythonJob_outputs")
    wg.add_task(
        add,
        name="add1",
        x=1,
        y=2,
        computer="localhost",
        command_info={"label": python_executable_path},
    )
    wg.add_task(
        decorted_multiply,
        name="multiply1",
        x=wg.tasks.add1.outputs.sum,
        y=3,
        computer="localhost",
        command_info={"label": python_executable_path},
    )
    # wg.submit(wait=True)
    wg.run()
    assert wg.tasks.add1.outputs.sum.value.value == 3
    assert wg.tasks.add1.outputs["diff"].value.value == -1
    assert wg.tasks.multiply1.outputs.result.value.value == 9
    # process_label and label
    assert wg.tasks.add1.node.process_label == "PythonJob<add1>"
    assert wg.tasks.add1.node.label == "add1"


def test_PythonJob_kwargs(fixture_localhost, python_executable_path):
    """Test function with kwargs."""

    @task.pythonjob()
    def add(x, y=1, **kwargs):
        x += y
        for value in kwargs.values():
            x += value
        return x

    wg = WorkGraph("test_PythonJob")
    wg.add_task(add, name="add1")
    wg.run(
        inputs={
            "add1": {
                "x": 1,
                "y": 2,
                "kwargs": {"m": 2, "n": 3},
                "computer": "localhost",
                "command_info": {"label": python_executable_path},
            },
        },
    )
    # data inside the kwargs should be serialized separately
    wg.process.inputs.tasks.add1.kwargs.m == 2
    assert wg.tasks.add1.outputs.result.value.value == 8
    # load the workgraph
    wg = WorkGraph.load(wg.pk)
    assert wg.tasks.add1.inputs["kwargs"]._value == {"m": 2, "n": 3}


def test_dynamic_inputs(fixture_localhost, python_executable_path) -> None:
    """Test dynamic inputs.
    For dynamic inputs, we allow the user to define the inputs manually.
    """

    @task.pythonjob()
    def add(**kwargs):
        return sum(kwargs.values())

    wg = WorkGraph("test_dynamic_inputs")
    wg.add_task(
        add,
        name="add1",
        x=np.array([1, 2]),
        y=np.array([3, 4]),
        command_info={"label": python_executable_path},
    )
    wg.run()
    assert (wg.tasks.add1.outputs.result.value.get_array() == np.array([4, 6])).all()


def test_PythonJob_namespace_output_input(fixture_localhost, python_executable_path):
    """Test function with namespace output and input."""

    # output namespace
    out = namespace(
        add_multiply=namespace(add=namespace(sum=Any, total=Any), multiply=Any),
        minus=Any,
    )

    @task.pythonjob
    def myfunc(x, y) -> out:
        return {
            "add_multiply": {"add": {"sum": x + y, "total": x + y}, "multiply": x * y},
            "minus": x - y,
        }

    # input namespace
    @task.pythonjob
    def myfunc2(x, y):
        add = x["add"]["sum"]
        multiply = x["multiply"]
        return y + add + multiply

    @task.pythonjob
    def myfunc3(x, y):
        return x["sum"] + y

    wg = WorkGraph("test_namespace_outputs")
    wg.add_task(myfunc, name="myfunc")
    wg.add_task(
        myfunc2,
        name="myfunc2",
        x=wg.tasks.myfunc.outputs.add_multiply,
    )
    wg.add_task(
        myfunc3,
        name="myfunc3",
        x=wg.tasks.myfunc.outputs.add_multiply.add,
    )

    inputs = {
        "myfunc": {
            "x": 1.0,
            "y": 2.0,
            "computer": "localhost",
            "command_info": {"label": python_executable_path},
        },
        "myfunc2": {
            "y": 3.0,
            "computer": "localhost",
            "command_info": {"label": python_executable_path},
        },
        "myfunc3": {
            "y": 4.0,
            "computer": "localhost",
            "command_info": {"label": python_executable_path},
        },
    }
    wg.run(inputs=inputs)
    assert wg.tasks.myfunc.outputs.add_multiply.add.sum.value == 3
    assert wg.tasks.myfunc.outputs.add_multiply.multiply.value == 2
    assert wg.tasks.myfunc2.outputs.result.value == 8
    assert wg.tasks.myfunc3.outputs.result.value == 7


def test_PythonJob_copy_files(fixture_localhost, python_executable_path):
    """Test function with copy files."""

    # define add task
    @task.pythonjob()
    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    # define multiply task
    @task.pythonjob()
    def multiply(x_folder_name, y_folder_name):
        with open(f"{x_folder_name}/result.txt", "r") as f:
            x = int(f.read())
        with open(f"{y_folder_name}/result.txt", "r") as f:
            y = int(f.read())
        return x * y

    wg = WorkGraph("test_PythonJob_parent_folder")
    wg.add_task(add, name="add1")
    wg.add_task(add, name="add2")
    wg.add_task(
        multiply,
        name="multiply",
    )
    wg.add_link(
        wg.tasks.add1.outputs.remote_folder,
        wg.tasks["multiply"].inputs["copy_files"],
    )
    wg.add_link(
        wg.tasks.add2.outputs.remote_folder,
        wg.tasks["multiply"].inputs["copy_files"],
    )
    wg.run(
        inputs={
            "add1": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
                "command_info": {"label": python_executable_path},
            },
            "add2": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
                "command_info": {"label": python_executable_path},
            },
            "multiply": {
                "x_folder_name": "add1_remote_folder",
                "y_folder_name": "add2_remote_folder",
                "computer": "localhost",
                "command_info": {"label": python_executable_path},
            },
        },
    )
    assert wg.tasks["multiply"].outputs.result.value.value == 25


def test_load_pythonjob(fixture_localhost, python_executable_path):
    """Test function with typing."""

    @task.pythonjob()
    def add(x: str, y: str) -> str:
        return x + y

    wg = WorkGraph("test_PythonJob")
    wg.add_task(add, name="add")

    wg.run(
        inputs={
            "add": {
                "x": "Hello, ",
                "y": "World!",
                "computer": "localhost",
                "command_info": {"label": python_executable_path},
            },
        },
        # wait=True,
    )
    assert wg.tasks.add.outputs.result.value.value == "Hello, World!"
    wg = WorkGraph.load(wg.pk)
    wg.tasks.add.inputs.x.value = "Hello, "
    wg.tasks.add.inputs["y"].value = "World!"


def test_exit_code(fixture_localhost, python_executable_path):
    """Test function with exit code."""

    def handle_negative_sum(task: Task):
        """Handle the failure code 410 of the `add`.
        Simply make the inputs positive by taking the absolute value.
        """
        # because the error_handler is ran by the engine
        # all the inputs are serialized to AiiDA data
        # therefore, we need use the `value` attribute to get the raw data
        task.set_inputs(
            {
                "x": abs(task.inputs.x.value.get_array()),
                "y": abs(task.inputs.y.value.get_array()),
            }
        )

        return "Run error handler: handle_negative_sum."

    @task.pythonjob(
        outputs=["sum"],
        error_handlers={
            "handle_negative_sum": {
                "executor": handle_negative_sum,
                "exit_codes": [410],
                "max_retries": 5,
            }
        },
    )
    def add(x: np.array, y: np.array) -> np.array:
        sum = x + y
        if (sum < 0).any():
            exit_code = {"status": 410, "message": "Some elements are negative"}
            return {"sum": sum, "exit_code": exit_code}
        return {"sum": sum}

    wg = WorkGraph("test_PythonJob")
    wg.add_task(
        add,
        name="add1",
        x=np.array([1, 1]),
        y=np.array([1, -2]),
        computer="localhost",
        command_info={"label": python_executable_path},
    )
    wg.run()
    # the first task should have exit status 410
    assert wg.process.base.links.get_outgoing().all()[0].node.exit_status == 410
    assert (
        wg.process.base.links.get_outgoing().all()[0].node.exit_message
        == "Some elements are negative"
    )
    # the final task should have exit status 0
    assert wg.tasks.add1.process.exit_status == 0
    assert (wg.tasks.add1.outputs.sum.value.get_array() == np.array([2, 3])).all()
