import aiida
from aiida_workgraph import WorkGraph, node
from typing import Any

aiida.load_profile()


def test_pythonjob_kwargs():
    """Test function with kwargs."""
    from aiida_workgraph import node, WorkGraph

    # define add node
    @node()
    def add(x, **kwargs):
        for value in kwargs.values():
            x += value
        return x

    wg = WorkGraph("test_pythonjob")
    wg.nodes.new(add, name="add", run_remotely=True)
    wg.run(
        inputs={
            "add": {
                "x": 1,
                "kwargs": {"y": 2, "z": 3},
                "computer": "localhost",
            },
        },
    )
    assert wg.nodes["add"].outputs["result"].value.value == 6


def test_pythonjob_typing():
    """Test function with typing."""
    from aiida_workgraph import node, WorkGraph
    from numpy import array

    # define add node
    @node()
    def add(x: array, y: array) -> array:
        return x + y

    # define multiply node
    @node()
    def multiply(x: Any, y: Any) -> Any:
        return x * y

    wg = WorkGraph("test_pythonjob")
    wg.nodes.new(add, name="add", run_remotely=True)
    wg.nodes.new(
        multiply, name="multiply", run_remotely=True, x=wg.nodes["add"].outputs[0]
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
    assert (wg.nodes["multiply"].outputs["result"].value.value == array([12, 20])).all()


def test_pythonjob_outputs():
    """Test function with multiple outputs."""

    @node(outputs=[["General", "sum"], ["General", "diff"]])
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    wg = WorkGraph("test_pythonjob_outputs")
    wg.nodes.new(
        add,
        name="add",
        x=1,
        y=2,
        run_remotely=True,
        #  code=code,
        computer="localhost",
    )
    wg.submit(wait=True)
    assert wg.nodes["add"].outputs["sum"].value.value == 3
    assert wg.nodes["add"].outputs["diff"].value.value == -1


def test_pythonjob_parent_folder():
    """Test function with parent folder."""
    from aiida_workgraph import WorkGraph, node
    from aiida import load_profile

    load_profile()

    # define add node
    @node()
    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    # define multiply node
    @node()
    def multiply(x, y):
        with open("parent_folder/result.txt", "r") as f:
            z = int(f.read())
        return x * y + z

    wg = WorkGraph("first_workflow")
    wg.nodes.new(add, name="add", run_remotely=True)
    wg.nodes.new(
        multiply,
        name="multiply",
        parent_folder=wg.nodes["add"].outputs["remote_folder"],
        run_remotely=True,
    )

    # ------------------------- Submit the calculation -------------------
    wg.submit(
        inputs={
            "add": {
                "x": 2,
                "y": 3,
                # "code": code,
                "computer": "localhost",
            },
            "multiply": {
                "x": 3,
                "y": 4,
                #  "code": code,
                "computer": "localhost",
            },
        },
        wait=True,
    )
    assert wg.nodes["multiply"].outputs["result"].value.value == 17


def test_pythonjob_upload_files():
    """Test function with upload files."""
    from aiida_workgraph import WorkGraph, node

    # create a temporary file "input.txt" in the current directory
    with open("input.txt", "w") as f:
        f.write("2")

    # create a temporary folder "inputs_folder" in the current directory
    # and add a file "another_input.txt" in the folder
    import os

    os.makedirs("inputs_folder", exist_ok=True)
    with open("inputs_folder/another_input.txt", "w") as f:
        f.write("3")

    # define add node
    @node()
    def add():
        with open("input.txt", "r") as f:
            a = int(f.read())
        with open("inputs_folder/another_input.txt", "r") as f:
            b = int(f.read())
        return a + b

    wg = WorkGraph("first_workflow")
    wg.nodes.new(add, name="add", run_remotely=True)

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
    assert wg.nodes["add"].outputs["result"].value.value == 5
