import aiida
from aiida_workgraph import WorkGraph, task
from typing import Any

aiida.load_profile()


def test_PythonTask_kwargs():
    """Test function with kwargs."""
    from aiida_workgraph import task, WorkGraph

    # define add task
    @task()
    def add(x, **kwargs):
        for value in kwargs.values():
            x += value
        return x

    wg = WorkGraph("test_PythonTask")
    wg.tasks.new(add, name="add", run_remotely=True)
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


def test_PythonTask_typing():
    """Test function with typing."""
    from aiida_workgraph import task, WorkGraph
    from numpy import array

    # define add task
    @task()
    def add(x: array, y: array) -> array:
        return x + y

    # define multiply task
    @task()
    def multiply(x: Any, y: Any) -> Any:
        return x * y

    wg = WorkGraph("test_PythonTask")
    wg.tasks.new(add, name="add", run_remotely=True)
    wg.tasks.new(
        multiply, name="multiply", run_remotely=True, x=wg.tasks["add"].outputs[0]
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


def test_PythonTask_outputs():
    """Test function with multiple outputs."""

    @task(outputs=[["General", "sum"], ["General", "diff"]])
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    wg = WorkGraph("test_PythonTask_outputs")
    wg.tasks.new(
        add,
        name="add",
        x=1,
        y=2,
        run_remotely=True,
        #  code=code,
        computer="localhost",
    )
    wg.submit(wait=True)
    assert wg.tasks["add"].outputs["sum"].value.value == 3
    assert wg.tasks["add"].outputs["diff"].value.value == -1


def test_PythonTask_parent_folder():
    """Test function with parent folder."""
    from aiida_workgraph import WorkGraph, task
    from aiida import load_profile

    load_profile()

    # define add task
    @task()
    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    # define multiply task
    @task()
    def multiply(x, y):
        with open("parent_folder/result.txt", "r") as f:
            z = int(f.read())
        return x * y + z

    wg = WorkGraph("test_PythonTask_parent_folder")
    wg.tasks.new(add, name="add", run_remotely=True)
    wg.tasks.new(
        multiply,
        name="multiply",
        parent_folder=wg.tasks["add"].outputs["remote_folder"],
        run_remotely=True,
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


def test_PythonTask_upload_files():
    """Test function with upload files."""
    from aiida_workgraph import WorkGraph, task

    # create a temporary file "input.txt" in the current directory
    with open("input.txt", "w") as f:
        f.write("2")

    # create a temporary folder "inputs_folder" in the current directory
    # and add a file "another_input.txt" in the folder
    import os

    os.makedirs("inputs_folder", exist_ok=True)
    with open("inputs_folder/another_input.txt", "w") as f:
        f.write("3")

    # define add task
    @task()
    def add():
        with open("input.txt", "r") as f:
            a = int(f.read())
        with open("inputs_folder/another_input.txt", "r") as f:
            b = int(f.read())
        return a + b

    wg = WorkGraph("test_PythonTask_upload_files")
    wg.tasks.new(add, name="add", run_remotely=True)

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


def test_PythonTask_copy_files():
    """Test function with copy files."""
    from aiida_workgraph import WorkGraph, task

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

    wg = WorkGraph("test_PythonTask_parent_folder")
    wg.tasks.new(add, name="add1", run_remotely=True)
    wg.tasks.new(add, name="add2", run_remotely=True)
    wg.tasks.new(
        multiply,
        name="multiply",
        run_remotely=True,
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
