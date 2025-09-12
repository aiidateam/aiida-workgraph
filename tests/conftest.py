import pytest
from aiida_workgraph import task, WorkGraph
from aiida.orm import Int, StructureData
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from typing import Callable, Any, Union
from aiida.orm import WorkflowNode
import time
import os
from aiida_workgraph.socket_spec import namespace

pytest_plugins = [
    "aiida.tools.pytest_fixtures",
]


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
    from aiida.orm import InstalledCode, load_code
    from aiida.common import NotExistent

    try:
        code = load_code("add@localhost")
    except NotExistent:
        code = InstalledCode(
            label="add",
            computer=fixture_localhost,
            filepath_executable="/bin/bash",
            default_calc_job_plugin="arithmetic.add",
        )
        code.store()
    return code


@pytest.fixture(scope="session")
def python_executable_path():
    return os.environ.get("PYTEST_PYTHONJOB_PYTHON_EXEC_PATH", "python3")


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
def wg_task() -> WorkGraph:
    """A workgraph with calcfunction."""

    wg = WorkGraph(name="test_debug_math")
    sumdiff1 = wg.add_task("workgraph.test_sum_diff", "sumdiff1", x=2, y=3)
    sumdiff2 = wg.add_task("workgraph.test_sum_diff", "sumdiff2", x=4)
    wg.add_link(sumdiff1.outputs[0], sumdiff2.inputs[1])
    return wg


@pytest.fixture
def wg_calcjob(add_code) -> WorkGraph:
    """A workgraph with calcjob."""

    print("add_code", add_code)

    wg = WorkGraph(name="test_debug_math")
    add1 = wg.add_task(ArithmeticAddCalculation, "add1", x=2, y=3, code=add_code)
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=4, code=add_code)
    wg.add_link(add1.outputs.sum, add2.inputs["y"])
    return wg


@pytest.fixture
def wg_workchain(add_code) -> WorkGraph:
    """A workgraph with workchain."""
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    wg = WorkGraph(name="test_debug_math")
    int1 = wg.add_task("workgraph.load_node", "int1", pk=Int(2).store().pk)
    int2 = wg.add_task("workgraph.load_node", "int2", pk=Int(3).store().pk)
    code1 = wg.add_task("workgraph.load_code", "code1", pk=add_code.pk)
    multiply_add1 = wg.add_task(MultiplyAddWorkChain, "multiply_add1", x=Int(4).store())
    multiply_add2 = wg.add_task(
        MultiplyAddWorkChain,
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

    @task.calcfunction()
    def add(
        x: Union[int, float], y: Union[int, float], t: Union[int, float] = 1.0
    ) -> Union[int, float]:

        time.sleep(t.value)
        return x + y

    return add


@pytest.fixture
def decorated_multiply() -> Callable:
    """Generate a decorated node for test."""

    @task.calcfunction
    def multiply(
        x: Union[int, float], y: Union[int, float], t: Union[int, float] = 1.0
    ) -> Union[int, float]:

        time.sleep(t.value)
        return x * y

    return multiply


@pytest.fixture
def decorated_sqrt() -> Callable:
    """Generate a decorated node for test."""

    @task.calcfunction
    def mysqrt(x: Union[int, float]) -> Union[int, float]:
        from math import sqrt

        return sqrt(x)

    return mysqrt


@pytest.fixture
def decorated_smaller_than() -> Callable:
    """Generate a decorated node for test."""

    # define smaller_than task
    @task()
    def smaller_than(x, y):
        return x < y

    return smaller_than


@pytest.fixture
def decorated_add_multiply(decorated_add, decorated_multiply) -> Callable:
    """Generate a decorated node for test."""

    @task.workfunction
    def add_multiply(x, y, z):
        """Add two numbers and multiply it with a third."""
        # we need use the calcfunction to get the result, instead of the wrapped function
        addition = decorated_add._func(x, y)
        product = decorated_multiply._func(addition, z)
        return {"result": product}

    return add_multiply


@pytest.fixture
def decorated_add_multiply_group(decorated_add, decorated_multiply) -> Callable:
    """Generate a decorated node for test."""

    @task.graph()
    def add_multiply_group(x, y, z):
        outputs1 = decorated_add(x=x, y=y)
        outputs2 = decorated_multiply(x=z, y=outputs1.result)
        return outputs2.result

    return add_multiply_group


@pytest.fixture
def decorated_namespace_sum_diff() -> Callable:
    """Generate a decorated node for test."""

    out = namespace(sum=Any, diff=Any, nested=namespace(diff=Any, sum=Any))

    @task
    def sum_diff(x, y, nested: namespace(x=Any, y=Any)) -> out:
        """Add two numbers and return the result."""
        return {
            "sum": x + y,
            "diff": x - y,
            "nested": {
                "diff": nested["x"] - nested["y"],
                "sum": nested["x"] + nested["y"],
            },
        }

    return sum_diff


@pytest.fixture
def structure_si() -> StructureData:
    from ase.build import bulk

    si = bulk("Si")
    structure_si = StructureData(ase=si)
    return structure_si


@pytest.fixture
def wg_engine(decorated_add, add_code) -> WorkGraph:
    """Use to test the engine."""
    code = add_code
    wg = WorkGraph(name="test_run_order")
    add0 = wg.add_task(ArithmeticAddCalculation, "add0", x=2, y=0, code=code)
    add1 = wg.add_task(decorated_add, "add1", x=2, y=1)
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=2, y=2, code=code)
    add3 = wg.add_task(decorated_add, "add3", x=2, y=3)
    add4 = wg.add_task(ArithmeticAddCalculation, "add4", x=2, y=4, code=code)
    add5 = wg.add_task(decorated_add, "add5", x=2, y=5)
    wg.add_link(add0.outputs.sum, add2.inputs.x)
    wg.add_link(add1.outputs[0], add3.inputs.x)
    wg.add_link(add3.outputs[0], add4.inputs.x)
    wg.add_link(add2.outputs.sum, add5.inputs.x)
    wg.add_link(add4.outputs.sum, add5.inputs["y"])
    return wg


@pytest.fixture
def create_process_node():
    """Return a process node."""

    def process_node(state="finished", exit_status=0):
        """Return a finished process node."""

        node = WorkflowNode()
        node.set_process_state(state)
        node.set_exit_status(exit_status)
        node.seal()
        node.store()
        return node

    return process_node


@pytest.fixture
def create_workgraph_process_node():
    """Return a process node."""

    def process_node(state="finished", exit_status=0):
        """Return a finished process node."""
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        process = WorkGraphEngine(
            inputs={
                "workgraph_data": {
                    "name": "test",
                    "state": "",
                    "tasks": {},
                    "links": [],
                }
            }
        )
        node = process.node
        node.set_process_state(state)
        node.set_exit_status(exit_status)
        node.seal()
        node.store()
        return node

    return process_node
