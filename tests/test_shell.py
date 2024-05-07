import aiida
from aiida_workgraph import WorkGraph
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData

aiida.load_profile()


def test_shell_date():
    # Create a code on the local computer
    cat_code = prepare_code("cat")
    # Create a workgraph
    wg = WorkGraph(name="test_shell_cat_with_file_arguments")
    job1 = wg.nodes.new(
        "AiiDAShell",
        code=cat_code,
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
    from aiida_shell.launch import prepare_code
    from aiida.orm import Int
    from aiida_shell.data import PickledData
    import os

    echo_code = prepare_code("echo")
    bc_code = prepare_code("bc")

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
        "AiiDAShell",
        name="job1",
        code=echo_code,
        arguments=["{x}", "+", "{y}"],
        nodes={
            "x": Int(2),
            "y": Int(3),
        },
    )
    # bc command to calculate the expression
    job2 = wg.nodes.new(
        "AiiDAShell",
        name="job2",
        code=bc_code,
        arguments=["{expression}"],
        filenames={"expression": "input.txt"},
        parser=PickledData(parser),
    )
    # add a "result" output socket from the parser
    job2.outputs.new("General", "result")
    # echo result + y expression
    job3 = wg.nodes.new(
        "AiiDAShell",
        name="job3",
        code=echo_code,
        arguments=["{result}", "*", "{z}"],
    )
    # bc command to calculate the expression
    job4 = wg.nodes.new(
        "AiiDAShell",
        name="job4",
        code=bc_code,
        arguments=["{expression}"],
        filenames={"expression": "input.txt"},
        parser=PickledData(parser),
    )
    # add a "result" output socket from the parser
    job4.outputs.new("General", "result")
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
