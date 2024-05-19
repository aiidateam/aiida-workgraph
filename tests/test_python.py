import aiida
from aiida_workgraph import WorkGraph, node

aiida.load_profile()


def test_python_job():
    """Test a simple python node."""
    from aiida_workgraph import node, WorkGraph

    code = aiida.orm.load_code("python@localhost")

    # define add node
    @node()
    def add(x, y):
        return x + y

    # define multiply node
    @node()
    def multiply(x, y):
        return x * y

    wg = WorkGraph("test_python_job")
    wg.nodes.new(add, name="add", on_remote=True)
    wg.nodes.new(
        multiply, name="multiply", on_remote=True, x=wg.nodes["add"].outputs[0]
    )
    wg.submit(
        inputs={
            "add": {"x": 2, "y": 3, "_code": code},
            "multiply": {"y": 4, "_code": code},
        },
        wait=True,
    )
    assert wg.nodes["multiply"].outputs["result"].value.value == 20


def test_python_job_outputs():
    """Test a simple python node."""
    code = aiida.orm.load_code("python@localhost")

    @node(outputs=[["General", "sum"], ["General", "diff"]])
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    wg = WorkGraph("test_python_job_outputs")
    wg.nodes.new(add, name="add", on_remote=True, x=1, y=2, _code=code)
    wg.submit(wait=True)
    assert wg.nodes["add"].outputs["sum"].value.value == 3
    assert wg.nodes["add"].outputs["diff"].value.value == -1
