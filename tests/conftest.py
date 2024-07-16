import pytest
from aiida_workgraph import task, WorkGraph
from aiida.engine import calcfunction, workfunction
from aiida.orm import Float, Int, StructureData
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from typing import Callable, Any, Union
import time

pytest_plugins = "aiida.tools.pytest_fixtures"


@pytest.fixture(scope="session", autouse=True)
def aiida_profile(aiida_config, aiida_profile_factory):
    """Create and load a profile with RabbitMQ as broker."""
    with aiida_profile_factory(aiida_config, broker_backend="core.rabbitmq") as profile:
        yield profile


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def add_code(fixture_localhost):
    from aiida.orm import InstalledCode

    code = InstalledCode(computer=fixture_localhost, filepath_executable="/bin/bash")
    code.store()
    return code


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return an ``InstalledCode`` instance configured to run calculations of given entry point on localhost."""

    def _fixture_code(entry_point_name):
        from aiida.common import exceptions
        from aiida.orm import InstalledCode, load_code

        label = f"test.{entry_point_name}"

        try:
            return load_code(label=label)
        except exceptions.NotExistent:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable="/bin/true",
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def wg_calcfunction() -> WorkGraph:
    """A workgraph with calcfunction."""

    wg = WorkGraph(name="test_debug_math")
    float1 = wg.add_task("AiiDANode", "float1", pk=Float(3.0).store().pk)
    sumdiff1 = wg.add_task("AiiDASumDiff", "sumdiff1", x=2)
    sumdiff2 = wg.add_task("AiiDASumDiff", "sumdiff2", x=4)
    sumdiff3 = wg.add_task("AiiDASumDiff", "sumdiff3", x=6)
    wg.add_link(float1.outputs[0], sumdiff1.inputs[1])
    wg.add_link(sumdiff1.outputs[0], sumdiff2.inputs[1])
    wg.add_link(sumdiff2.outputs[0], sumdiff3.inputs[1])
    return wg


@pytest.fixture
def wg_calcjob(add_code) -> WorkGraph:
    """A workgraph with calcjob."""

    print("add_code", add_code)

    wg = WorkGraph(name="test_debug_math")
    int1 = wg.add_task("AiiDANode", "int1", pk=Int(3).store().pk)
    code1 = wg.add_task("AiiDACode", "code1", pk=add_code.pk)
    add1 = wg.add_task(ArithmeticAddCalculation, "add1", x=Int(2).store())
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=Int(4).store())
    add3 = wg.add_task(ArithmeticAddCalculation, "add3", x=Int(4).store())
    wg.add_link(code1.outputs[0], add1.inputs["code"])
    wg.add_link(int1.outputs[0], add1.inputs["y"])
    wg.add_link(code1.outputs[0], add2.inputs["code"])
    wg.add_link(add1.outputs["sum"], add2.inputs["y"])
    wg.add_link(code1.outputs[0], add3.inputs["code"])
    wg.add_link(add2.outputs["sum"], add3.inputs["y"])
    return wg


@pytest.fixture
def wg_workchain(add_code) -> WorkGraph:
    """A workgraph with workchain."""

    wg = WorkGraph(name="test_debug_math")
    int1 = wg.add_task("AiiDANode", "int1", pk=Int(2).store().pk)
    int2 = wg.add_task("AiiDANode", "int2", pk=Int(3).store().pk)
    code1 = wg.add_task("AiiDACode", "code1", pk=add_code.pk)
    multiply_add1 = wg.add_task(
        "AiiDAArithmeticMultiplyAdd", "multiply_add1", x=Int(4).store()
    )
    multiply_add2 = wg.add_task(
        "AiiDAArithmeticMultiplyAdd",
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    wg.add_link(code1.outputs[0], multiply_add1.inputs["code"])
    wg.add_link(int1.outputs[0], multiply_add1.inputs["y"])
    wg.add_link(int2.outputs[0], multiply_add1.inputs["z"])
    wg.add_link(code1.outputs[0], multiply_add2.inputs["code"])
    wg.add_link(multiply_add1.outputs[0], multiply_add2.inputs["z"])
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

    # define compare task
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

    @task.graph_builder(outputs=[{"name": "result", "from": "multiply1.result"}])
    def add_multiply_group(x, y, z):
        wg = WorkGraph("add_multiply_group")
        add1 = wg.add_task(decorated_add, name="add1", x=x, y=y)
        multiply = wg.add_task(decorated_multiply, name="multiply1", x=z)
        # link the output of a task to the input of another task
        wg.add_link(add1.outputs[0], multiply.inputs["y"])
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
    structure1 = wg.add_task("AiiDAStructure", "structure1")
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
def wg_engine(decorated_add, add_code) -> WorkGraph:
    """Use to test the engine."""
    code = add_code
    x = Int(2)
    wg = WorkGraph(name="test_run_order")
    add0 = wg.add_task(ArithmeticAddCalculation, "add0", x=x, y=Int(0), code=code)
    add0.set({"metadata.options.sleep": 15})
    add1 = wg.add_task(decorated_add, "add1", x=x, y=Int(1), t=Int(1))
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=x, y=Int(2), code=code)
    add2.set({"metadata.options.sleep": 1})
    add3 = wg.add_task(decorated_add, "add3", x=x, y=Int(3), t=Int(1))
    add4 = wg.add_task(ArithmeticAddCalculation, "add4", x=x, y=Int(4), code=code)
    add4.set({"metadata.options.sleep": 1})
    add5 = wg.add_task(decorated_add, "add5", x=x, y=Int(5), t=Int(1))
    wg.add_link(add0.outputs["sum"], add2.inputs["x"])
    wg.add_link(add1.outputs[0], add3.inputs["x"])
    wg.add_link(add3.outputs[0], add4.inputs["x"])
    wg.add_link(add2.outputs["sum"], add5.inputs["x"])
    wg.add_link(add4.outputs["sum"], add5.inputs["y"])
    return wg
