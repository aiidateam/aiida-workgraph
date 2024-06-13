import pytest
from aiida_workgraph import task, WorkGraph
from aiida.engine import calcfunction, workfunction
from aiida.orm import Float, Int, load_code, StructureData
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from typing import Callable, Any, Union
import time


@pytest.fixture
def wg_calcfunction() -> WorkGraph:
    """A workgraph with calcfunction."""

    wg = WorkGraph(name="test_debug_math")
    float1 = wg.tasks.new("AiiDANode", "float1", pk=Float(3.0).store().pk)
    sumdiff1 = wg.tasks.new("AiiDASumDiff", "sumdiff1", x=2)
    sumdiff2 = wg.tasks.new("AiiDASumDiff", "sumdiff2", x=4)
    sumdiff3 = wg.tasks.new("AiiDASumDiff", "sumdiff3", x=6)
    wg.links.new(float1.outputs[0], sumdiff1.inputs[1])
    wg.links.new(sumdiff1.outputs[0], sumdiff2.inputs[1])
    wg.links.new(sumdiff2.outputs[0], sumdiff3.inputs[1])
    return wg


@pytest.fixture
def wg_calcjob() -> WorkGraph:
    """A workgraph with calcjob."""

    code = load_code("add@localhost")
    wg = WorkGraph(name="test_debug_math")
    int1 = wg.tasks.new("AiiDANode", "int1", pk=Int(3).store().pk)
    code1 = wg.tasks.new("AiiDACode", "code1", pk=code.pk)
    add1 = wg.tasks.new(ArithmeticAddCalculation, "add1", x=Int(2).store())
    add2 = wg.tasks.new(ArithmeticAddCalculation, "add2", x=Int(4).store())
    add3 = wg.tasks.new(ArithmeticAddCalculation, "add3", x=Int(4).store())
    wg.links.new(code1.outputs[0], add1.inputs["code"])
    wg.links.new(int1.outputs[0], add1.inputs["y"])
    wg.links.new(code1.outputs[0], add2.inputs["code"])
    wg.links.new(add1.outputs["sum"], add2.inputs["y"])
    wg.links.new(code1.outputs[0], add3.inputs["code"])
    wg.links.new(add2.outputs["sum"], add3.inputs["y"])
    return wg


@pytest.fixture
def wg_workchain() -> WorkGraph:
    """A workgraph with workchain."""

    code = load_code("add@localhost")
    wg = WorkGraph(name="test_debug_math")
    int1 = wg.tasks.new("AiiDANode", "int1", pk=Int(2).store().pk)
    int2 = wg.tasks.new("AiiDANode", "int2", pk=Int(3).store().pk)
    code1 = wg.tasks.new("AiiDACode", "code1", pk=code.pk)
    multiply_add1 = wg.tasks.new(
        "AiiDAArithmeticMultiplyAdd", "multiply_add1", x=Int(4).store()
    )
    multiply_add2 = wg.tasks.new(
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
def decorated_normal_add() -> Callable:
    """Generate a decorated node for test."""

    @task()
    def add(x: Any, y: Any) -> Any:
        return x + y

    return add


@pytest.fixture
def decorated_add() -> Callable:
    """Generate a decorated node for test."""

    @calcfunction
    def add(
        x: Union[int, float], y: Union[int, float], t: Union[int, float] = 1.0
    ) -> Union[int, float]:

        time.sleep(t.value)
        return x + y

    return add


@pytest.fixture
def decorated_multiply() -> Callable:
    """Generate a decorated node for test."""

    @calcfunction
    def multiply(
        x: Union[int, float], y: Union[int, float], t: Union[int, float] = 1.0
    ) -> Union[int, float]:

        time.sleep(t.value)
        return x * y

    return multiply


@pytest.fixture
def decorated_sqrt() -> Callable:
    """Generate a decorated node for test."""

    @calcfunction
    def mysqrt(x: Union[int, float]) -> Union[int, float]:
        from math import sqrt

        return sqrt(x)

    return mysqrt


@pytest.fixture
def decorated_compare() -> Callable:
    """Generate a decorated node for test."""

    # define compare node
    @task()
    def compare(x, y):
        return x < y

    return compare


@pytest.fixture
def decorated_add_multiply(decorated_add, decorated_multiply) -> Callable:
    """Generate a decorated node for test."""

    @workfunction
    def add_multiply(x, y, z):
        """Add two numbers and multiply it with a third."""
        addition = decorated_add(x, y)
        product = decorated_multiply(addition, z)
        return {"result": product}

    return add_multiply


@pytest.fixture
def decorated_add_multiply_group(decorated_add, decorated_multiply) -> Callable:
    """Generate a decorated node for test."""

    @task.graph_builder(outputs=[["multiply1.result", "result"]])
    def add_multiply_group(x, y, z):
        wg = WorkGraph("add_multiply_group")
        add1 = wg.tasks.new(decorated_add, name="add1", x=x, y=y)
        multiply = wg.tasks.new(decorated_multiply, name="multiply1", x=z)
        # link the output of int node to the input of add node
        wg.links.new(add1.outputs[0], multiply.inputs["y"])
        return wg

    return add_multiply_group


@pytest.fixture
def structure_si() -> StructureData:
    from ase.build import bulk

    si = bulk("Si")
    structure_si = StructureData(ase=si)
    return structure_si


@pytest.fixture
def wg_structure_si() -> WorkGraph:
    wg = WorkGraph(name="test_structure")
    structure1 = wg.tasks.new("AiiDAStructure", "structure1")
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
def wg_engine(decorated_add) -> WorkGraph:
    """Use to test the engine."""
    code = load_code("add@localhost")
    x = Int(2)
    wg = WorkGraph(name="test_run_order")
    add0 = wg.tasks.new(ArithmeticAddCalculation, "add0", x=x, y=Int(0), code=code)
    add0.set({"metadata.options.sleep": 15})
    add1 = wg.tasks.new(decorated_add, "add1", x=x, y=Int(1), t=Int(1))
    add2 = wg.tasks.new(ArithmeticAddCalculation, "add2", x=x, y=Int(2), code=code)
    add2.set({"metadata.options.sleep": 1})
    add3 = wg.tasks.new(decorated_add, "add3", x=x, y=Int(3), t=Int(1))
    add4 = wg.tasks.new(ArithmeticAddCalculation, "add4", x=x, y=Int(4), code=code)
    add4.set({"metadata.options.sleep": 1})
    add5 = wg.tasks.new(decorated_add, "add5", x=x, y=Int(5), t=Int(1))
    wg.links.new(add0.outputs["sum"], add2.inputs["x"])
    wg.links.new(add1.outputs[0], add3.inputs["x"])
    wg.links.new(add3.outputs[0], add4.inputs["x"])
    wg.links.new(add2.outputs["sum"], add5.inputs["x"])
    wg.links.new(add4.outputs["sum"], add5.inputs["y"])
    return wg
