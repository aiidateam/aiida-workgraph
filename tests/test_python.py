import pytest
from aiida_workgraph import WorkGraph, task
from typing import Any


@pytest.mark.usefixtures("started_daemon_client")
def test_PythonJob_kwargs(fixture_localhost):
    """Test function with kwargs."""

    def add(x, **kwargs):
        for value in kwargs.values():
            x += value
        return x

    wg = WorkGraph("test_PythonJob")
    wg.tasks.new("PythonJob", function=add, name="add")
    wg.run(
        inputs={
            "add": {
                "x": 1,
                "kwargs": {"y": 2, "z": 3},
                "computer": "localhost",
            },
        },
    )
    assert wg.tasks["add"].outputs["result"].value.value == 6


def test_PythonJob_typing(fixture_localhost):
    """Test function with typing."""
    from numpy import array

    def add(x: array, y: array) -> array:
        return x + y

    def multiply(x: Any, y: Any) -> Any:
        return x * y

    wg = WorkGraph("test_PythonJob")
    wg.tasks.new("PythonJob", function=add, name="add")
    wg.tasks.new(
        "PythonJob", function=multiply, name="multiply", x=wg.tasks["add"].outputs[0]
    )
    #
    metadata = {
        "options": {
            "custom_scheduler_commands": "# test",
            # "custom_scheduler_commands": 'module load anaconda\nconda activate py3.11\n',
        }
    }
    wg.run(
        inputs={
            "add": {
                "x": array([1, 2]),
                "y": array([2, 3]),
                "computer": "localhost",
                "metadata": metadata,
            },
            "multiply": {"y": 4, "computer": "localhost", "metadata": metadata},
        },
        # wait=True,
    )
    assert (wg.tasks["multiply"].outputs["result"].value.value == array([12, 20])).all()


def test_PythonJob_outputs(fixture_localhost):
    """Test function with multiple outputs."""

    @task(
        outputs=[
            {"identifier": "Any", "name": "sum"},
            {"identifier": "Any", "name": "diff"},
        ]
    )
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    wg = WorkGraph("test_PythonJob_outputs")
    wg.tasks.new(
        "PythonJob",
        function=add,
        name="add",
        x=1,
        y=2,
        #  code=code,
        computer="localhost",
    )
    # wg.submit(wait=True)
    wg.run()
    assert wg.tasks["add"].outputs["sum"].value.value == 3
    assert wg.tasks["add"].outputs["diff"].value.value == -1


def test_PythonJob_namespace_output(fixture_localhost):
    """Test function with namespace output and input."""

    # output namespace
    @task(
        outputs=[
            {
                "name": "add_multiply",
                "identifier": "Namespace",
            },
            {
                "name": "add_multiply.add",
                "identifier": "Namespace",
            },
            {"name": "minus"},
        ]
    )
    def myfunc(x, y):
        add = {"order1": x + y, "order2": x * x + y * y}
        return {
            "add_multiply": {"add": add, "multiply": x * y},
            "minus": x - y,
        }

    wg = WorkGraph("test_namespace_outputs")
    wg.tasks.new("PythonJob", function=myfunc, name="myfunc")
    wg.submit(
        wait=True,
        inputs={
            "myfunc": {
                "x": 1.0,
                "y": 2.0,
                "computer": "localhost",
            }
        },
    )
    assert wg.tasks["myfunc"].outputs["add_multiply"].value.add.order1.value == 3
    assert wg.tasks["myfunc"].outputs["add_multiply"].value.add.order2.value == 5
    assert wg.tasks["myfunc"].outputs["add_multiply"].value.multiply.value == 2


def test_PythonJob_namespace_output_input(fixture_localhost):
    """Test function with namespace output and input."""

    # output namespace
    @task(
        outputs=[
            {"identifier": "Namespace", "name": "add_multiply"},
            {"name": "add_multiply.add"},
            {"name": "add_multiply.multiply"},
            {"name": "minus"},
        ]
    )
    def myfunc(x, y):
        return {
            "add_multiply": {"add": x + y, "multiply": x * y},
            "minus": x - y,
        }

    # input namespace
    @task()
    def myfunc2(x, y):
        add = x["add"]
        multiply = x["multiply"]
        return y + add + multiply

    @task()
    def myfunc3(x, y):
        return x + y

    wg = WorkGraph("test_namespace_outputs")
    wg.tasks.new("PythonJob", function=myfunc, name="myfunc")
    wg.tasks.new(
        "PythonJob",
        function=myfunc2,
        name="myfunc2",
        x=wg.tasks["myfunc"].outputs["add_multiply"],
    )
    wg.tasks.new(
        "PythonJob",
        function=myfunc3,
        name="myfunc3",
        x=wg.tasks["myfunc"].outputs["add_multiply.add"],
    )

    inputs = {
        "myfunc": {
            "x": 1.0,
            "y": 2.0,
            "computer": "localhost",
        },
        "myfunc2": {
            "y": 3.0,
            "computer": "localhost",
        },
        "myfunc3": {
            "y": 4.0,
            "computer": "localhost",
        },
    }
    wg.run(inputs=inputs)
    assert wg.tasks["myfunc"].outputs["add_multiply"].value.add.value == 3
    assert wg.tasks["myfunc"].outputs["add_multiply"].value.multiply.value == 2
    assert wg.tasks["myfunc2"].outputs["result"].value.value == 8
    assert wg.tasks["myfunc3"].outputs["result"].value.value == 7


def test_PythonJob_parent_folder(fixture_localhost):
    """Test function with parent folder."""

    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    def multiply(x, y):
        with open("parent_folder/result.txt", "r") as f:
            z = int(f.read())
        return x * y + z

    wg = WorkGraph("test_PythonJob_parent_folder")
    wg.tasks.new("PythonJob", function=add, name="add")
    wg.tasks.new(
        "PythonJob",
        function=multiply,
        name="multiply",
        parent_folder=wg.tasks["add"].outputs["remote_folder"],
    )

    # ------------------------- Submit the calculation -------------------
    wg.submit(
        inputs={
            "add": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
            },
            "multiply": {
                "x": 3,
                "y": 4,
                "computer": "localhost",
            },
        },
        wait=True,
    )
    assert wg.tasks["multiply"].outputs["result"].value.value == 17


def test_PythonJob_upload_files(fixture_localhost):
    """Test function with upload files."""

    # create a temporary file "input.txt" in the current directory
    with open("input.txt", "w") as f:
        f.write("2")

    # create a temporary folder "inputs_folder" in the current directory
    # and add a file "another_input.txt" in the folder
    import os

    os.makedirs("inputs_folder", exist_ok=True)
    with open("inputs_folder/another_input.txt", "w") as f:
        f.write("3")

    def add():
        with open("input.txt", "r") as f:
            a = int(f.read())
        with open("inputs_folder/another_input.txt", "r") as f:
            b = int(f.read())
        return a + b

    wg = WorkGraph("test_PythonJob_upload_files")
    wg.tasks.new("PythonJob", function=add, name="add")

    # ------------------------- Submit the calculation -------------------
    # we need use full path to the file
    input_file = os.path.abspath("input.txt")
    input_folder = os.path.abspath("inputs_folder")

    wg.run(
        inputs={
            "add": {
                "computer": "localhost",
                "upload_files": {
                    "input.txt": input_file,
                    "inputs_folder": input_folder,
                },
            },
        },
    )
    # wait=True)
    assert wg.tasks["add"].outputs["result"].value.value == 5


def test_PythonJob_copy_files(fixture_localhost):
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
    wg.tasks.new("PythonJob", function=add, name="add1")
    wg.tasks.new("PythonJob", function=add, name="add2")
    wg.tasks.new(
        "PythonJob",
        function=multiply,
        name="multiply",
    )
    wg.links.new(
        wg.tasks["add1"].outputs["remote_folder"],
        wg.tasks["multiply"].inputs["copy_files"],
    )
    wg.links.new(
        wg.tasks["add2"].outputs["remote_folder"],
        wg.tasks["multiply"].inputs["copy_files"],
    )
    # ------------------------- Submit the calculation -------------------
    wg.submit(
        inputs={
            "add1": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
            },
            "add2": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
            },
            "multiply": {
                "x_folder_name": "add1_remote_folder",
                "y_folder_name": "add2_remote_folder",
                "computer": "localhost",
            },
        },
        wait=True,
    )
    assert wg.tasks["multiply"].outputs["result"].value.value == 25


def test_PythonJob_retrieve_files(fixture_localhost):
    """Test retrieve files."""

    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    wg = WorkGraph("test_PythonJob_retrieve_files")
    wg.tasks.new("PythonJob", function=add, name="add")
    # ------------------------- Submit the calculation -------------------
    wg.submit(
        inputs={
            "add": {
                "x": 2,
                "y": 3,
                "computer": "localhost",
                "metadata": {
                    "options": {
                        "additional_retrieve_list": ["result.txt"],
                    }
                },
            },
        },
        wait=True,
    )
    assert (
        "result.txt" in wg.tasks["add"].outputs["retrieved"].value.list_object_names()
    )


def test_data_serializer(fixture_localhost):
    from ase import Atoms
    from ase.build import bulk

    @task()
    def make_supercell(atoms: Atoms, dim: int) -> Atoms:
        """Scale the structure by the given scales."""
        return atoms * (dim, dim, dim)

    atoms = bulk("Si")

    wg = WorkGraph("test_PythonJob_retrieve_files")
    wg.tasks.new(
        "PythonJob", function=make_supercell, atoms=atoms, dim=2, name="make_supercell"
    )
    # ------------------------- Submit the calculation -------------------
    wg.submit(wait=True)
    assert (
        wg.tasks["make_supercell"].outputs["result"].value.value.get_chemical_formula()
        == "Si16"
    )
