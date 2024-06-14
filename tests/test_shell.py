import aiida
from aiida_workgraph import WorkGraph
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData

aiida.load_profile()


def test_shell_command():
    """Test the ShellTask with command as a string."""
    wg = WorkGraph(name="test_shell_command")
    job1 = wg.tasks.new(
        "ShellTask",
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
    """Test the ShellTask with code."""
    cat_code = prepare_code("cat")
    wg = WorkGraph(name="test_shell_code")
    job1 = wg.tasks.new(
        "ShellTask",
        command=cat_code,
        arguments=["{file_a}", "{file_b}"],
        nodes={
            "file_a": SinglefileData.from_string("string a"),
            "file_b": SinglefileData.from_string("string b"),
        },
    )
    wg.submit(wait=True)
    assert job1.node.outputs.stdout.get_content() == "string astring b"


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
        "ShellTask",
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
        "ShellTask",
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
        "ShellTask",
        name="job3",
        command="echo",
        arguments=["{result}", "*", "{z}"],
        nodes={"result": job2.outputs["result"], "z": Int(4)},
    )
    # bc command to calculate the expression
    job4 = wg.tasks.new(
        "ShellTask",
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
