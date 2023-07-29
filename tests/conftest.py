import pytest
from aiida_worktree import node, WorkTree
from aiida.engine import calcfunction, workfunction
from aiida.orm import Float, Int, Bool


@pytest.fixture
def nt_calcfunction():
    """A worktree with calcfunction."""

    nt = WorkTree(name="test_debug_math")
    float1 = nt.nodes.new("AiiDANode", "float1", value=Float(3.0).store())
    sumdiff1 = nt.nodes.new("AiiDASumDiff", "sumdiff1", x=2)
    sumdiff2 = nt.nodes.new("AiiDASumDiff", "sumdiff2", x=4)
    sumdiff3 = nt.nodes.new("AiiDASumDiff", "sumdiff3", x=6)
    nt.links.new(float1.outputs[0], sumdiff1.inputs[1])
    nt.links.new(sumdiff1.outputs[0], sumdiff2.inputs[1])
    nt.links.new(sumdiff2.outputs[0], sumdiff3.inputs[1])
    return nt


@pytest.fixture
def nt_calcjob():
    """A worktree with calcjob."""
    from aiida.orm import load_code

    code = load_code("add@localhost")
    nt = WorkTree(name="test_debug_math")
    int1 = nt.nodes.new("AiiDANode", "int1", value=Int(3).store())
    code1 = nt.nodes.new("AiiDACode", "code1", value=code.pk)
    add1 = nt.nodes.new("AiiDAArithmeticAdd", "add1", x=Int(2).store())
    add2 = nt.nodes.new("AiiDAArithmeticAdd", "add2", x=Int(4).store())
    nt.links.new(code1.outputs[0], add1.inputs[0])
    nt.links.new(int1.outputs[0], add1.inputs[2])
    nt.links.new(code1.outputs[0], add2.inputs[0])
    nt.links.new(add1.outputs[0], add2.inputs[2])
    return nt


@pytest.fixture
def nt_workchain():
    """A worktree with workchain."""
    from aiida.orm import load_code

    code = load_code("add@localhost")
    nt = WorkTree(name="test_debug_math")
    int1 = nt.nodes.new("AiiDANode", "int1", value=Int(2).store())
    int2 = nt.nodes.new("AiiDANode", "int2", value=Int(3).store())
    code1 = nt.nodes.new("AiiDACode", "code1", value=code.pk)
    multiply_add1 = nt.nodes.new(
        "AiiDAArithmeticMultiplyAdd", "multiply_add1", x=Int(4).store()
    )
    multiply_add2 = nt.nodes.new(
        "AiiDAArithmeticMultiplyAdd",
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    nt.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    nt.links.new(int1.outputs[0], multiply_add1.inputs["y"])
    nt.links.new(int2.outputs[0], multiply_add1.inputs["z"])
    nt.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    nt.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    return nt


@pytest.fixture
def decorated_add():
    """Generate a decorated node for test."""

    @node()
    @calcfunction
    def add(x, y):
        return x + y

    return add


@pytest.fixture
def decorated_multiply():
    """Generate a decorated node for test."""

    @node()
    @calcfunction
    def multiply(x, y):
        return x * y

    return multiply


@pytest.fixture
def decorated_sqrt():
    """Generate a decorated node for test."""

    @node()
    @calcfunction
    def mysqrt(x):
        from math import sqrt

        return sqrt(x)

    return mysqrt


@pytest.fixture
def decorated_compare():
    """Generate a decorated node for test."""

    # define compare node
    @node()
    @calcfunction
    def compare(x, y):
        return Bool(x < y)

    return compare


@pytest.fixture
def decorated_add_multiply(decorated_add, decorated_multiply):
    """Generate a decorated node for test."""

    @node()
    @workfunction
    def add_multiply(x, y, z):
        """Add two numbers and multiply it with a third."""
        addition = decorated_add(x, y)
        product = decorated_multiply(addition, z)
        return {"result": product}

    return add_multiply


@pytest.fixture
def decorated_add_multiply_group(decorated_add, decorated_multiply):
    """Generate a decorated node for test."""

    @node.group(outputs=[["multiply", "result", "result"]])
    def add_multiply_group(x, y, z):
        nt = WorkTree("add_multiply_group")
        add1 = nt.nodes.new(decorated_add, name="add1", x=x, y=y)
        multiply = nt.nodes.new(decorated_multiply, name="multiply", x=z)
        # link the output of int node to the input of add node
        nt.links.new(add1.outputs[0], multiply.inputs["y"])
        return nt

    return add_multiply_group


@pytest.fixture
def build_workchain():
    """Generate a decorated node for test."""

    from aiida_worktree import build_node

    ndata = {"path": "aiida.workflows.arithmetic.multiply_add.MultiplyAddWorkChain"}
    multiply_add = build_node(ndata)

    return multiply_add


@pytest.fixture
def structure_si():

    from aiida.orm import StructureData
    from ase.build import bulk

    si = bulk("Si")
    structure_si = StructureData(ase=si)
    return structure_si


@pytest.fixture
def nt_structure_si():

    nt = WorkTree(name="test_structure")
    structure1 = nt.nodes.new("AiiDAStructure", "structure1")
    data = {
        "cell": [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
        "kinds": [{"mass": 28.085, "name": "Si", "symbols": ["Si"], "weights": [1.0]}],
        "pbc1": True,
        "pbc2": True,
        "pbc3": True,
        "sites": [
            {"kind_name": "Si", "position": [0.0, 0.0, 0.0]},
            {"kind_name": "Si", "position": [1.3575, 1.3575, 1.3575]},
        ],
    }
    structure1.set(data)
    return nt
