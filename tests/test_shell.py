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
