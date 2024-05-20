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
    wg.nodes.new("PythonJob", function=add, name="add")
    wg.nodes.new(
        "PythonJob", function=multiply, name="multiply", x=wg.nodes["add"].outputs[0]
    )
    #
    metadata = {
        "options": {
            "custom_scheduler_commands": "# test",
        }
    }
    wg.submit(
        inputs={
            "add": {"x": 2, "y": 3, "code": code, "metadata": metadata},
            "multiply": {"y": 4, "code": code, "metadata": metadata},
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
    wg.nodes.new("PythonJob", function=add, name="add", x=1, y=2, code=code)
    wg.submit(wait=True)
    assert wg.nodes["add"].outputs["sum"].value.value == 3
    assert wg.nodes["add"].outputs["diff"].value.value == -1


def test_python_job_parent_folder():
    from aiida_workgraph import WorkGraph, node
    from aiida import orm, load_profile

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
    wg.nodes.new("PythonJob", function=add, name="add")
    wg.nodes.new(
        "PythonJob",
        function=multiply,
        name="multiply",
        parent_folder=wg.nodes["add"].outputs["remote_folder"],
    )

    # ------------------------- Submit the calculation -------------------
    code = orm.load_code("python@localhost")
    wg.submit(
        inputs={
            "add": {"x": 2, "y": 3, "code": code},
            "multiply": {"x": 3, "y": 4, "code": code},
        },
        wait=True,
    )
    assert wg.nodes["multiply"].outputs["result"].value.value == 17
