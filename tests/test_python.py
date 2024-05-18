import aiida
from aiida_workgraph import WorkGraph, node

aiida.load_profile()


def test_python():
    """Test a simple python node."""
    code = aiida.orm.load_code("python@localhost")

    @node()
    def add(x, y):
        return x + y

    wg = WorkGraph("test_python_node")
    wg.nodes.new(add, name="add", x=1, y=2, _code=code)
    wg.submit(wait=True)
    assert wg.nodes["add"].outputs["result"].value == 3
