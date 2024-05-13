from aiida_workgraph import worknode, WorkGraph
from aiida import load_profile, orm
from typing import Callable

load_profile()


def test_for(decorated_add: Callable, decorated_multiply: Callable) -> None:
    # Create a WorkGraph will loop the a sequence
    @worknode.graph_builder(outputs=[["context.total", "result"]])
    def add_multiply_for(sequence):
        wg = WorkGraph("add_multiply_for")
        # tell the engine that this is a `for` workgraph
        wg.workgraph_type = "FOR"
        # the sequence to be iter
        wg.sequence = sequence
        # set a context variable before running.
        wg.context = {"total": 0}
        multiply1 = wg.nodes.new(
            decorated_multiply, name="multiply1", x="{{ i }}", y=orm.Int(2)
        )
        add1 = wg.nodes.new(decorated_add, name="add1", x="{{ total }}")
        # update the context variable
        add1.to_context = [["result", "total"]]
        wg.links.new(multiply1.outputs["result"], add1.inputs["y"])
        # don't forget to return the workgraph
        return wg

    # -----------------------------------------
    wg = WorkGraph("test_for")
    for1 = wg.nodes.new(add_multiply_for, sequence=range(5))
    add1 = wg.nodes.new(decorated_add, name="add1", y=orm.Int(1))
    wg.links.new(for1.outputs["result"], add1.inputs["x"])
    wg.submit(wait=True, timeout=200)
    assert add1.node.outputs.result.value == 21
