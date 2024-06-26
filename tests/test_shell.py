import aiida
from aiida_workgraph import WorkGraph
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData

aiida.load_profile()


def test_shell_command():
    """Test the ShellJob with command as a string."""
    wg = WorkGraph(name="test_shell_command")
    job1 = wg.tasks.new(
        "ShellJob",
        command="cat",
        resolve_command=True,
        arguments=["{file_a}", "{file_b}"],
        nodes={
            "file_a": SinglefileData.from_string("string a"),
            "file_b": SinglefileData.from_string("string b"),
        },
    )
    wg.submit(wait=True)
    assert job1.node.outputs.stdout.get_content() == "string astring b"


def test_shell_code():
    """Test the ShellJob with code."""
    cat_code = prepare_code("cat")
    wg = WorkGraph(name="test_shell_code")
    job1 = wg.tasks.new(
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


def test_shell_set():
    """Set the nodes during/after the creation of the task."""
    wg = WorkGraph(name="test_shell_set")
    echo_task = wg.tasks.new(
        "ShellJob",
        name="echo",
        command="cp",
        arguments=["{file}", "copied_file"],
        nodes={"file": SinglefileData.from_string("1 5 1")},
        outputs=["copied_file"],
    )

    cat_task = wg.tasks.new(
        "ShellJob",
        name="cat",
        command="cat",
        arguments=["{input}"],
        nodes={"input": None},
    )
    wg.links.new(echo_task.outputs["copied_file"], cat_task.inputs["nodes.input"])
    wg.submit(wait=True)
    assert cat_task.outputs["stdout"].value.get_content() == "1 5 1"


def test_shell_workflow():
    from aiida_workgraph import WorkGraph
    from aiida.orm import Int
    from aiida_shell.data import PickledData

    def parser(self, dirpath):
        from aiida.orm import Int

        return {"result": Int((dirpath / "stdout").read_text().strip())}

    # Create a workgraph
    wg = WorkGraph(name="shell_add_mutiply_workflow")
    # echo x + y expression
    job1 = wg.tasks.new(
        "ShellJob",
        name="job1",
        command="echo",
        arguments=["{x}", "+", "{y}"],
        nodes={
            "x": Int(2),
            "y": Int(3),
        },
    )
    # bc command to calculate the expression
    job2 = wg.tasks.new(
        "ShellJob",
        name="job2",
        command="bc",
        arguments=["{expression}"],
        nodes={"expression": job1.outputs["stdout"]},
        parser=PickledData(parser),
        parser_outputs=[
            ["General", "result"]
        ],  # add a "result" output socket from the parser
    )
    # echo result + y expression
    job3 = wg.tasks.new(
        "ShellJob",
        name="job3",
        command="echo",
        arguments=["{result}", "*", "{z}"],
        nodes={"result": job2.outputs["result"], "z": Int(4)},
    )
    # bc command to calculate the expression
    job4 = wg.tasks.new(
        "ShellJob",
        name="job4",
        command="bc",
        arguments=["{expression}"],
        nodes={"expression": job3.outputs["stdout"]},
        parser=PickledData(parser),
        parser_outputs=[
            ["General", "result"]
        ],  # add a "result" output socket from the parser
    )
    # there is a bug in aiida-shell, the following line will raise an error
    # https://github.com/sphuber/aiida-shell/issues/91
    # wg.submit(wait=True, timeout=200)
    wg.run()
    assert job4.outputs["result"].value.value == 20
