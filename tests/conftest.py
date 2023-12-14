import pytest
from aiida_worktree import node, WorkTree, build_node
from aiida.orm import Float, Int, Bool, load_code


@pytest.fixture
def arithmetic_add():
    """Generate a node for test."""

    arithmetic_add = build_node(
        {"path": "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"}
    )
    return arithmetic_add


@pytest.fixture
def wt_calcfunction():
    """A worktree with calcfunction."""

    wt = WorkTree(name="test_debug_math")
    float1 = wt.nodes.new("AiiDANode", "float1", value=Float(3.0).store())
    sumdiff1 = wt.nodes.new("AiiDASumDiff", "sumdiff1", x=2)
    sumdiff2 = wt.nodes.new("AiiDASumDiff", "sumdiff2", x=4)
    sumdiff3 = wt.nodes.new("AiiDASumDiff", "sumdiff3", x=6)
    wt.links.new(float1.outputs[0], sumdiff1.inputs[1])
    wt.links.new(sumdiff1.outputs[0], sumdiff2.inputs[1])
    wt.links.new(sumdiff2.outputs[0], sumdiff3.inputs[1])
    return wt


@pytest.fixture
def wt_calcjob(arithmetic_add):
    """A worktree with calcjob."""

    code = load_code("add@localhost")
    wt = WorkTree(name="test_debug_math")
    int1 = wt.nodes.new("AiiDANode", "int1", value=Int(3).store())
    code1 = wt.nodes.new("AiiDACode", "code1", value=code.pk)
    add1 = wt.nodes.new(arithmetic_add, "add1", x=Int(2).store())
    add2 = wt.nodes.new(arithmetic_add, "add2", x=Int(4).store())
    add3 = wt.nodes.new(arithmetic_add, "add3", x=Int(4).store())
    wt.links.new(code1.outputs[0], add1.inputs["code"])
    wt.links.new(int1.outputs[0], add1.inputs["y"])
    wt.links.new(code1.outputs[0], add2.inputs["code"])
    wt.links.new(add1.outputs["sum"], add2.inputs["y"])
    wt.links.new(code1.outputs[0], add3.inputs["code"])
    wt.links.new(add2.outputs["sum"], add3.inputs["y"])
    return wt


@pytest.fixture
def wt_workchain():
    """A worktree with workchain."""

    code = load_code("add@localhost")
    wt = WorkTree(name="test_debug_math")
    int1 = wt.nodes.new("AiiDANode", "int1", value=Int(2).store())
    int2 = wt.nodes.new("AiiDANode", "int2", value=Int(3).store())
    code1 = wt.nodes.new("AiiDACode", "code1", value=code.pk)
    multiply_add1 = wt.nodes.new(
        "AiiDAArithmeticMultiplyAdd", "multiply_add1", x=Int(4).store()
    )
    multiply_add2 = wt.nodes.new(
        "AiiDAArithmeticMultiplyAdd",
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    wt.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    wt.links.new(int1.outputs[0], multiply_add1.inputs["y"])
    wt.links.new(int2.outputs[0], multiply_add1.inputs["z"])
    wt.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    wt.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    return wt


@pytest.fixture
def decorated_normal_add():
    """Generate a decorated node for test."""

    @node()
    def add(x, y):
        return x + y

    return add


@pytest.fixture
def decorated_add():
    """Generate a decorated node for test."""

    @node.calcfunction()
    def add(x, y, t=1):
        import time

        time.sleep(t.value)
        return x + y

    return add


@pytest.fixture
def decorated_multiply():
    """Generate a decorated node for test."""

    @node.calcfunction()
    def multiply(x, y, t=1):
        import time

        time.sleep(t.value)
        return x * y

    return multiply


@pytest.fixture
def decorated_sqrt():
    """Generate a decorated node for test."""

    @node.calcfunction()
    def mysqrt(x):
        from math import sqrt

        return sqrt(x)

    return mysqrt


@pytest.fixture
def decorated_compare():
    """Generate a decorated node for test."""

    # define compare node
    @node.calcfunction()
    def compare(x, y):
        return Bool(x < y)

    return compare


@pytest.fixture
def decorated_add_multiply(decorated_add, decorated_multiply):
    """Generate a decorated node for test."""

    @node.workfunction()
    def add_multiply(x, y, z):
        """Add two numbers and multiply it with a third."""
        addition = decorated_add(x, y)
        product = decorated_multiply(addition, z)
        return {"result": product}

    return add_multiply


@pytest.fixture
def decorated_add_multiply_group(decorated_add, decorated_multiply):
    """Generate a decorated node for test."""

    @node.group(outputs=[["multiply.result", "result"]])
    def add_multiply_group(x, y, z):
        wt = WorkTree("add_multiply_group")
        add1 = wt.nodes.new(decorated_add, name="add1", x=x, y=y)
        multiply = wt.nodes.new(decorated_multiply, name="multiply", x=z)
        # link the output of int node to the input of add node
        wt.links.new(add1.outputs[0], multiply.inputs["y"])
        return wt

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
def wt_structure_si():

    wt = WorkTree(name="test_structure")
    structure1 = wt.nodes.new("AiiDAStructure", "structure1")
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
    return wt


@pytest.fixture
def wt_engine(arithmetic_add):
    """Use to test the engine."""
    code = load_code("add@localhost")
    x = Int(2)
    wt = WorkTree(name="test_run_order")
    adds = []
    for i in range(6):
        temp = wt.nodes.new(arithmetic_add, f"add{i}", x=x, y=Int(i), code=code)
        if i == 0:
            temp.set({"metadata.options.sleep": 15})
        else:
            temp.set({"metadata.options.sleep": 1})
        adds.append(temp)
    wt.links.new(adds[0].outputs["sum"], adds[2].inputs["x"])
    wt.links.new(adds[1].outputs["sum"], adds[3].inputs["x"])
    wt.links.new(adds[3].outputs["sum"], adds[4].inputs["x"])
    wt.links.new(adds[2].outputs["sum"], adds[5].inputs["x"])
    wt.links.new(adds[4].outputs["sum"], adds[5].inputs["y"])
    return wt
