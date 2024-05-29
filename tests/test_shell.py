import aiida
from aiida_workgraph import WorkGraph
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData

aiida.load_profile()


def test_shell_command():
    """Test the ShellJob with command as a string."""
    wg = WorkGraph(name="test_shell_command")
    job1 = wg.nodes.new(
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
    job1 = wg.nodes.new(
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


def test_shell_workflow():
    from aiida_workgraph import node, WorkGraph
    from aiida.orm import Int
    from aiida_shell.data import PickledData
    import os

    def parser(self, dirpath):
        from aiida.orm import Int

        return {"result": Int((dirpath / "stdout").read_text().strip())}

    @node()
    def prepare_bc_nodes(file):
        """Prepare the nodes for the bc calculation."""
        return {"result": {"expression": file}}

    @node()
    def prepare_echo_nodes(result, z):
        """Prepare the nodes for the echo calculation."""
        return {
            "result": {
                "result": result,
                "z": z,
            }
        }

    # Create a workgraph
    wg = WorkGraph(name="shell_add_mutiply_workflow")
    # echo x + y expression
    job1 = wg.nodes.new(
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
    job2 = wg.nodes.new(
        "ShellJob",
        name="job2",
        command="bc",
        arguments=["{expression}"],
        filenames={"expression": "input.txt"},
        parser=PickledData(parser),
        add_outputs=[
            ["General", "result"]
        ],  # add a "result" output socket from the parser
    )
    # echo result + y expression
    job3 = wg.nodes.new(
        "ShellJob",
        name="job3",
        command="echo",
        arguments=["{result}", "*", "{z}"],
    )
    # bc command to calculate the expression
    job4 = wg.nodes.new(
        "ShellJob",
        name="job4",
        command="bc",
        arguments=["{expression}"],
        filenames={"expression": "input.txt"},
        parser=PickledData(parser),
        add_outputs=[
            ["General", "result"]
        ],  # add a "result" output socket from the parser
    )
    # prepare the nodes for the bc calculation
    nodes1 = wg.nodes.new(prepare_bc_nodes, name="prepare_bc_nodes1")
    nodes2 = wg.nodes.new(prepare_echo_nodes, name="prepare_echo_nodes", z=Int(4))
    nodes3 = wg.nodes.new(prepare_bc_nodes, name="prepare_bc_nodes2")

    wg.links.new(job1.outputs["stdout"], nodes1.inputs["file"])
    wg.links.new(nodes1.outputs[0], job2.inputs["nodes"])
    wg.links.new(job2.outputs["result"], nodes2.inputs["result"])
    wg.links.new(nodes2.outputs["result"], job3.inputs["nodes"])
    wg.links.new(job3.outputs["stdout"], nodes3.inputs["file"])
    wg.links.new(nodes3.outputs[0], job4.inputs["nodes"])
    # there is a bug in aiida-shell, the following line will raise an error
    # https://github.com/sphuber/aiida-shell/issues/91
    # wg.submit(wait=True, timeout=200)
    wg.run()
    print("state: ", wg.state)
    os.system(f"verdi process report {wg.pk}")
    assert job4.outputs["result"].value.value == 20
