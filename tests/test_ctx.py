from aiida_workgraph import WorkGraph, task
from typing import Callable
from aiida.orm import Float, ArrayData
from aiida_workgraph import socket_spec as spec
from typing import Any
import numpy as np
import pytest


def test_ctx_update():
    wg = WorkGraph()
    wg.ctx = {"x": 1, "data.x": 5}
    assert wg.ctx.x.value == 1
    wg.update_ctx({"data.y": 10})
    assert wg.ctx.data.y.value == 10
    assert wg.ctx.x._link_limit > 1


def test_workgraph_ctx(decorated_add: Callable) -> None:
    """Set/get data to/from context."""

    wg = WorkGraph(name="test_workgraph_ctx")
    # create a array data object to test if it can be set to context
    # the workgraph should be able to serialize it
    array = ArrayData()
    array.set_array("matrix", np.array([[1, 2], [3, 4]]))
    wg.ctx = {"x": Float(2), "data.y": Float(3), "array": array}
    add1 = wg.add_task(decorated_add, "add1", x=wg.ctx.x, y=wg.ctx.data.y)
    wg.add_task(
        "workgraph.set_context", name="set_ctx1", key="x", value=add1.outputs.result
    )
    get_ctx1 = wg.add_task("workgraph.get_context", name="get_ctx1", key="x")
    # test the task can wait for another task
    get_ctx1.waiting_on.add(add1)
    add2 = wg.add_task(decorated_add, "add2", x=get_ctx1.outputs.result, y=1)
    wg.run()
    assert add2.outputs.result.value == 6


@pytest.mark.usefixtures("started_daemon_client")
def test_task_update_ctx(decorated_add: Callable) -> None:
    """Set/get data to/from context."""

    wg = WorkGraph(name="test_node_set_ctx")
    add1 = wg.add_task(decorated_add, "add1", x=Float(2).store(), y=Float(3).store())
    with pytest.raises(
        AttributeError, match="TaskSocketNamespace: add1.outputs has no attribute"
    ):
        wg.update_ctx({"sum": add1.outputs.resul})
    wg.update_ctx({"sum": add1.outputs.result})
    add2 = wg.add_task(decorated_add, "add2", x=add1.outputs[0], y=wg.ctx.sum)
    wg.run()
    assert add2.outputs.result.value == 10


def test_task_update_nested_ctx():
    @task
    def add(x, y) -> spec.namespace(results=spec.namespace(sum=Any, product=Any)):
        return {"results": {"sum": x + y, "product": x * y}}

    wg = WorkGraph()
    wg.add_task(add, name="add", x=1, y=2)
    wg.update_ctx({"data.sum": wg.tasks.add.outputs.results.sum})
    wg.outputs.sum = wg.ctx.data.sum
    wg.outputs.add = {}
    wg.outputs.add.sum = wg.ctx.data.sum
    wg.run()
    assert wg.process.outputs.sum.value == 3
    assert wg.process.outputs.add.sum.value == 3
