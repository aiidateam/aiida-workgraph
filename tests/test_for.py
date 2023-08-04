from aiida_worktree import node, WorkTree
from aiida import load_profile, orm

load_profile()


def test_for(decorated_add, decorated_multiply):
    # Create a WorkTree will loop the a sequence
    @node.group(outputs=[["ctx", "total", "result"]])
    def add_multiply_for(sequence):
        wt = WorkTree("add_multiply_for")
        # tell the engine that this is a `for` worktree
        wt.is_for = True
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
        # don't forget to return the worktree
        return wt

    # -----------------------------------------
    wt = WorkTree("test_for")
    for1 = wt.nodes.new(add_multiply_for, sequence=range(5))
    add1 = wt.nodes.new(decorated_add, y=orm.Int(1))
    wt.links.new(for1.outputs[0], add1.inputs[0])
    wt.submit(wait=True)
    assert add1.node.outputs.result.value == 21
