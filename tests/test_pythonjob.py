import pytest
from aiida_workgraph import WorkGraph, task, Task, TaskPool
from typing import Any
import numpy as np


def test_to_dict():
    """Test the to_dict method.
    PythonJobTask override the `to_dict` method to serialize the raw
    data to AiiDA data."""
    from ase import Atoms
    from ase.build import bulk
    from aiida import orm

    @task.pythonjob()
    def make_supercell(atoms: Atoms, dim: int = 2) -> Atoms:
        """Scale the structure by the given scales."""
        return atoms * (dim, dim, dim)

    atoms = bulk("Si")
    structure = orm.StructureData(ase=atoms)
    wg = WorkGraph("test_PythonJob_retrieve_files")
    # atoms will be converted to AtomsData automatically
    wg.add_task(
        make_supercell,
        atoms=atoms,
        dim=2,
        name="make_supercell_1",
    )
    # structure will be converted to AtomsData automatically
    wg.add_task(
        make_supercell,
        atoms=structure,
        dim=2,
        name="make_supercell_2",
    )
    data = wg.tasks.make_supercell_1.to_dict()
    assert isinstance(
        data["inputs"]["sockets"]["atoms"]["property"]["value"],
        orm.Data,
    )
    data = wg.tasks.make_supercell_1.to_dict()
    assert isinstance(
        data["inputs"]["sockets"]["atoms"]["property"]["value"],
        orm.Data,
    )


@pytest.mark.usefixtures("started_daemon_client")
def test_decorator(fixture_localhost, python_executable_path):
    """Test decorator."""

    @task.pythonjob(
        outputs=[
            {"identifier": "workgraph.any", "name": "sum"},
            {"identifier": "workgraph.any", "name": "diff"},
        ]
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

    def add(x, y=1, **kwargs):
        x += y
        for value in kwargs.values():
            x += value
        return x

    wg = WorkGraph("test_PythonJob")
    wg.add_task(TaskPool.workgraph.pythonjob, function=add, name="add1")
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
    wg.process.inputs.workgraph_data.tasks.add1.inputs.sockets.kwargs.sockets.m.property.value == 2
    assert wg.tasks.add1.outputs.result.value.value == 8
    # load the workgraph
    wg = WorkGraph.load(wg.pk)
    assert wg.tasks.add1.inputs["kwargs"]._value == {"m": 2, "n": 3}


def test_dynamic_inputs() -> None:
    """Test dynamic inputs.
    For dynamic inputs, we allow the user to define the inputs manually.
    """

    @task.pythonjob()
    def add(**kwargs):
        return sum(kwargs.values())

    wg = WorkGraph("test_dynamic_inputs")
    wg.add_task(add, name="add1", x=np.array([1, 2]), y=np.array([3, 4]))
    wg.run()
    assert (wg.tasks.add1.outputs.result.value.value == np.array([4, 6])).all()


def test_PythonJob_namespace_output_input(fixture_localhost, python_executable_path):
    """Test function with namespace output and input."""

    # output namespace
    @task(
        outputs=[
            {"name": "add_multiply", "identifier": "workgraph.namespace"},
            {
                "name": "add_multiply.add",
                "identifier": "workgraph.namespace",
                "metadata": {"dynamic": True},
            },
            {"name": "add_multiply.multiply"},
            {"name": "minus"},
        ]
    )
    def myfunc(x, y):
        return {
            "add_multiply": {"add": {"sum": x + y, "total": x + y}, "multiply": x * y},
            "minus": x - y,
        }

    # input namespace
    @task()
    def myfunc2(x, y):
        add = x["add"]["sum"]
        multiply = x["multiply"]
        return y + add + multiply

    @task()
    def myfunc3(x, y):
        return x["sum"] + y

    wg = WorkGraph("test_namespace_outputs")
    wg.add_task(TaskPool.workgraph.pythonjob, function=myfunc, name="myfunc")
    wg.add_task(
        TaskPool.workgraph.pythonjob,
        function=myfunc2,
        name="myfunc2",
        x=wg.tasks.myfunc.outputs.add_multiply,
    )
    wg.add_task(
        TaskPool.workgraph.pythonjob,
        function=myfunc3,
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
    @task()
    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    # define multiply task
    @task()
    def multiply(x_folder_name, y_folder_name):
        with open(f"{x_folder_name}/result.txt", "r") as f:
            x = int(f.read())
        with open(f"{y_folder_name}/result.txt", "r") as f:
            y = int(f.read())
        return x * y

    wg = WorkGraph("test_PythonJob_parent_folder")
    wg.add_task(TaskPool.workgraph.pythonjob, function=add, name="add1")
    wg.add_task(TaskPool.workgraph.pythonjob, function=add, name="add2")
    wg.add_task(
        TaskPool.workgraph.pythonjob,
        function=multiply,
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
    # ------------------------- Submit the calculation -------------------
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

        task.set(
            {
                "x": abs(task.inputs.x.value),
                "y": abs(task.inputs["y"].value),
            }
        )

        return "Run error handler: handle_negative_sum."

    @task.pythonjob(
        outputs=[{"name": "sum"}],
        error_handlers=[
            {"handler": handle_negative_sum, "exit_codes": [410], "max_retries": 5}
        ],
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
    assert wg.tasks.add1.node.exit_status == 0
    assert (wg.tasks.add1.outputs.sum.value.value == np.array([2, 3])).all()
