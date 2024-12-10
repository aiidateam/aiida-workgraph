import pytest
from aiida_workgraph import WorkGraph, task
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData, load_computer, Int


def test_prepare_for_shell_task_nonexistent():
    """Check that the `ValueError` raised by `aiida-shell` for a non-extistent executable is captured by WorkGraph."""
    from aiida_workgraph.engine.utils import prepare_for_shell_task

    task = {"name": "test"}
    inputs = {"command": "abc42"}
    with pytest.raises(ValueError, match="failed to determine the absolute path"):
        prepare_for_shell_task(task=task, inputs=inputs)


@pytest.mark.usefixtures("started_daemon_client")
def test_shell_command(fixture_localhost):
    """Test the ShellJob with command as a string."""
    wg = WorkGraph(name="test_shell_command")
    job1 = wg.add_task(
        "ShellJob",
        command="cat",
        resolve_command=True,
        arguments=["{file_a}", "{file_b}"],
        nodes={
            "file_a": SinglefileData.from_string("string a"),
            "file_b": SinglefileData.from_string("string b"),
        },
    )
    # also check if we can set the computer explicitly
    job1.set({"metadata.computer": load_computer("localhost")})
    wg.submit(wait=True)
    assert job1.node.outputs.stdout.get_content() == "string astring b"


def test_shell_code():
    """Test the ShellJob with code."""
    cat_code = prepare_code("cat")
    wg = WorkGraph(name="test_shell_code")
    job1 = wg.add_task(
        "ShellJob",
        command=cat_code,
        arguments=["{file_a}", "{file_b}"],
        nodes={
            "file_a": SinglefileData.from_string("string a"),
            "file_b": SinglefileData.from_string("string b"),
        },
    )
    wg.submit(wait=True)
    assert job1.node.outputs.stdout.get_content() == "string astring b"


def test_dynamic_port():
    """Set the nodes during/after the creation of the task."""
    wg = WorkGraph(name="test_dynamic_port")
    echo_task = wg.add_task(
        "ShellJob",
        name="echo",
        command="cp",
        arguments=["{file}", "copied_file"],
        nodes={"file": SinglefileData.from_string("1 5 1")},
        outputs=["copied_file"],
    )

    cat_task = wg.add_task(
        "ShellJob",
        name="cat",
        command="cat",
        arguments=["{input}"],
        nodes={"input1": None, "input2": Int(2), "input3": echo_task.outputs["_wait"]},
    )
    wg.add_link(echo_task.outputs["copied_file"], cat_task.inputs["nodes.input1"])
    # task will create input for each item in the dynamic port (nodes)
    assert "nodes.input1" in cat_task.inputs
    assert "nodes.input2" in cat_task.inputs
    # if the value of the item is a Socket, then it will create a link, and pop the item
    assert "nodes.input3" in cat_task.inputs
    assert cat_task.inputs["nodes"].socket_value == {"input2": Int(2)}


@pytest.mark.usefixtures("started_daemon_client")
def test_shell_graph_builder():
    """Test the ShellJob inside a graph builder.
    And the parser is also defined in the graph builder."""
    from aiida.orm import Int

    @task.graph_builder(outputs=[{"name": "result", "from": "job2.result"}])
    def add_multiply(x, y):
        """Add two numbers and multiply the result by 2."""
        # define the parser function
        def parser(dirpath):
            from aiida.orm import Int

            return {"result": Int((dirpath / "stdout").read_text().strip())}

        wg = WorkGraph()
        # echo x + y expression
        job1 = wg.add_task(
            "ShellJob",
            name="job1",
            command="echo",
            arguments=["{x}", "+", "{y}"],
            nodes={
                "x": x,
                "y": y,
            },
        )
        # bc command to calculate the expression
        wg.add_task(
            "ShellJob",
            name="job2",
            command="bc",
            arguments=["{expression}"],
            nodes={"expression": job1.outputs["stdout"]},
            parser=parser,
            parser_outputs=[
                {"identifier": "workgraph.any", "name": "result"}
            ],  # add a "result" output socket from the parser
        )
        return wg

    wg = WorkGraph(name="test_shell_graph_builder")
    add_multiply1 = wg.add_task(add_multiply, x=Int(2), y=Int(3))
    wg.submit(wait=True)
    assert add_multiply1.outputs["result"].socket_value.value == 5
