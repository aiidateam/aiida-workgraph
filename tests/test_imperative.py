from aiida_workgraph import task
from aiida_workgraph.engine.imperative.imperative import (
    WorkGraphImperativeEngine,
    wait_for,
)
from aiida.engine import run
import pytest


@task.pythonjob()
def add(x, y):
    return x + y


@task.pythonjob()
def multiply(x, y):
    return x * y


@pytest.mark.usefixtures("started_daemon_client")
def test_if(fixture_localhost):
    async def if_flow(x, y):
        add_result = add(x, y)
        await wait_for(add_result)
        add_result._node.graph.update()
        if add_result.result.value > 10:
            multiply_result = multiply(1, add_result.result)
        else:
            multiply_result = multiply(-1, add_result.result)
        return {"sum": add_result.result, "multiply": multiply_result.result}

    results = run(
        WorkGraphImperativeEngine,
        inputs={
            "workgraph_data": {
                "name": "if_flow",
                "flow": if_flow,
                "function_inputs": {"x": 3, "y": 4},
            }
        },
    )
    assert results["sum"].value == 7


@pytest.mark.usefixtures("started_daemon_client")
def test_while(fixture_localhost):
    async def add_multiply(x, y):
        outputs = add(x, y)
        await wait_for(outputs)
        while outputs.result.value < 10:
            outputs1 = add(outputs.result, 1)
            outputs = multiply(outputs1.result, y=2)
            await wait_for(outputs)
        return {"sum": outputs.result}

    results = run(
        WorkGraphImperativeEngine,
        inputs={
            "workgraph_data": {
                "name": "add_multiply",
                "flow": add_multiply,
                "function_inputs": {"x": 3, "y": 4},
            }
        },
    )
    assert results["sum"].value == 16
