from aiida_workgraph import node, WorkGraph
from aiida import load_profile, orm

load_profile()


def test_for(decorated_add, decorated_multiply):
    # Create a WorkGraph will loop the a sequence
    @node.group(outputs=[["ctx.total", "result"]])
    def add_multiply_for(sequence):
        wt = WorkGraph("add_multiply_for")
        # tell the engine that this is a `for` workgraph
        wt.workgraph_type = "FOR"
        # the sequence to be iter
        wt.sequence = sequence
        # set a context variable before running.
        wt.ctx = {"total": 0}
        multiply1 = wt.nodes.new(
            decorated_multiply, name="multiply1", x="{{ i }}", y=orm.Int(2)
        )
        add1 = wt.nodes.new(decorated_add, name="add1", x="{{ total }}")
        # update the context variable
        add1.to_ctx = [["result", "total"]]
        wt.links.new(multiply1.outputs[0], add1.inputs[1])
        # don't forget to return the workgraph
        return wt

    # -----------------------------------------
    wt = WorkGraph("test_for")
    for1 = wt.nodes.new(add_multiply_for, sequence=range(5))
    add1 = wt.nodes.new(decorated_add, y=orm.Int(1))
    wt.links.new(for1.outputs[0], add1.inputs[0])
    wt.submit(wait=True, timeout=200)
    assert add1.node.outputs.result.value == 21
