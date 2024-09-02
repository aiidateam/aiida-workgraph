from aiida_workgraph import WorkGraph
from typing import Callable
from aiida.orm import Float, ArrayData
import numpy as np


def test_workgraph_ctx(decorated_add: Callable) -> None:
    """Set/get data to/from context."""

    wg = WorkGraph(name="test_workgraph_ctx")
    # create a array data object to test if it can be set to context
    # the workgraph should be able to serialize it
    array = ArrayData()
    array.set_array("matrix", np.array([[1, 2], [3, 4]]))
    wg.context = {"x": Float(2), "data.y": Float(3), "array": array}
    add1 = wg.add_task(decorated_add, "add1", x="{{ x }}", y="{{ data.y }}")
    wg.add_task(
        "workgraph.to_context", name="to_ctx1", key="x", value=add1.outputs["result"]
    )
    from_ctx1 = wg.add_task("workgraph.from_context", name="from_ctx1", key="x")
    add2 = wg.add_task(decorated_add, "add2", x=from_ctx1.outputs["result"], y=1)
    wg.run()
    assert add2.outputs["result"].value == 6


def test_node_to_ctx(decorated_add: Callable) -> None:
    """Set/get data to/from context."""

    wg = WorkGraph(name="test_node_to_ctx")
    add1 = wg.add_task(decorated_add, "add1", x=Float(2).store(), y=Float(3).store())
    try:
        add1.set_context({"resul": "sum"})
    except ValueError as e:
        assert str(e) == "Keys {'resul'} are not in the outputs of this task."
    add1.set_context({"result": "sum"})
    add2 = wg.add_task(decorated_add, "add2", y="{{ sum }}")
    wg.add_link(add1.outputs[0], add2.inputs["x"])
    wg.submit(wait=True)
    assert add2.outputs["result"].value == 10
