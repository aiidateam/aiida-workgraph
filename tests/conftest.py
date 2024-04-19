import pytest
from aiida_workgraph import node, WorkGraph, build_node
from aiida.orm import Float, Int, load_code


@pytest.fixture
def arithmetic_add():
    """Generate a node for test."""

    arithmetic_add = build_node(
        "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"
    )
    return arithmetic_add


@pytest.fixture
def wg_calcfunction():
    """A workgraph with calcfunction."""

    wg = WorkGraph(name="test_debug_math")
    float1 = wg.nodes.new("AiiDANode", "float1", value=Float(3.0).store())
    sumdiff1 = wg.nodes.new("AiiDASumDiff", "sumdiff1", x=2)
    sumdiff2 = wg.nodes.new("AiiDASumDiff", "sumdiff2", x=4)
    sumdiff3 = wg.nodes.new("AiiDASumDiff", "sumdiff3", x=6)
    wg.links.new(float1.outputs[0], sumdiff1.inputs[1])
    wg.links.new(sumdiff1.outputs[0], sumdiff2.inputs[1])
    wg.links.new(sumdiff2.outputs[0], sumdiff3.inputs[1])
    return wg


@pytest.fixture
def wg_calcjob(arithmetic_add):
    """A workgraph with calcjob."""

    code = load_code("add@localhost")
    wg = WorkGraph(name="test_debug_math")
    int1 = wg.nodes.new("AiiDANode", "int1", value=Int(3).store())
    code1 = wg.nodes.new("AiiDACode", "code1", value=code.pk)
    add1 = wg.nodes.new(arithmetic_add, "add1", x=Int(2).store())
    add2 = wg.nodes.new(arithmetic_add, "add2", x=Int(4).store())
    add3 = wg.nodes.new(arithmetic_add, "add3", x=Int(4).store())
    wg.links.new(code1.outputs[0], add1.inputs["code"])
    wg.links.new(int1.outputs[0], add1.inputs["y"])
    wg.links.new(code1.outputs[0], add2.inputs["code"])
    wg.links.new(add1.outputs["sum"], add2.inputs["y"])
    wg.links.new(code1.outputs[0], add3.inputs["code"])
    wg.links.new(add2.outputs["sum"], add3.inputs["y"])
    return wg


@pytest.fixture
def wg_workchain():
    """A workgraph with workchain."""

    code = load_code("add@localhost")
    wg = WorkGraph(name="test_debug_math")
    int1 = wg.nodes.new("AiiDANode", "int1", value=Int(2).store())
    int2 = wg.nodes.new("AiiDANode", "int2", value=Int(3).store())
    code1 = wg.nodes.new("AiiDACode", "code1", value=code.pk)
    multiply_add1 = wg.nodes.new(
        "AiiDAArithmeticMultiplyAdd", "multiply_add1", x=Int(4).store()
    )
    multiply_add2 = wg.nodes.new(
        "AiiDAArithmeticMultiplyAdd",
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    wg.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    wg.links.new(int1.outputs[0], multiply_add1.inputs["y"])
    wg.links.new(int2.outputs[0], multiply_add1.inputs["z"])
    wg.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    wg.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    return wg


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
    @node()
    def compare(x, y):
        return x < y

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

    @node.group(outputs=[["multiply1.result", "result"]])
    def add_multiply_group(x, y, z):
        wg = WorkGraph("add_multiply_group")
        add1 = wg.nodes.new(decorated_add, name="add1", x=x, y=y)
        multiply = wg.nodes.new(decorated_multiply, name="multiply1", x=z)
        # link the output of int node to the input of add node
        wg.links.new(add1.outputs[0], multiply.inputs["y"])
        return wg

    return add_multiply_group


@pytest.fixture
def build_workchain():
    """Generate a decorated node for test."""

    from aiida_workgraph import build_node
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    multiply_add = build_node(MultiplyAddWorkChain)

    return multiply_add


@pytest.fixture
def structure_si():

    from aiida.orm import StructureData
    from ase.build import bulk

    si = bulk("Si")
    structure_si = StructureData(ase=si)
    return structure_si


@pytest.fixture
def wg_structure_si():

    wg = WorkGraph(name="test_structure")
    structure1 = wg.nodes.new("AiiDAStructure", "structure1")
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
    return wg


@pytest.fixture
def wg_engine(decorated_add, arithmetic_add):
    """Use to test the engine."""
    code = load_code("add@localhost")
    x = Int(2)
    wg = WorkGraph(name="test_run_order")
    add0 = wg.nodes.new(arithmetic_add, "add0", x=x, y=Int(0), code=code)
    add0.set({"metadata.options.sleep": 15})
    add1 = wg.nodes.new(decorated_add, "add1", x=x, y=Int(1), t=Int(1))
    add2 = wg.nodes.new(arithmetic_add, "add2", x=x, y=Int(2), code=code)
    add2.set({"metadata.options.sleep": 1})
    add3 = wg.nodes.new(decorated_add, "add3", x=x, y=Int(3), t=Int(1))
    add4 = wg.nodes.new(arithmetic_add, "add4", x=x, y=Int(4), code=code)
    add4.set({"metadata.options.sleep": 1})
    add5 = wg.nodes.new(decorated_add, "add5", x=x, y=Int(5), t=Int(1))
    wg.links.new(add0.outputs["sum"], add2.inputs["x"])
    wg.links.new(add1.outputs[0], add3.inputs["x"])
    wg.links.new(add3.outputs[0], add4.inputs["x"])
    wg.links.new(add2.outputs["sum"], add5.inputs["x"])
    wg.links.new(add4.outputs["sum"], add5.inputs["y"])
    return wg
