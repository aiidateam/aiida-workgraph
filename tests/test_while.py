from aiida_workgraph import node, WorkGraph
from aiida import load_profile, orm

load_profile()


def test_while(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    @node.group(outputs=[["ctx.n", "result"]])
    def my_while(n, limit):
        wt = WorkGraph("while_workgraph")
        wt.workgraph_type = "WHILE"
        wt.conditions = ["compare1.result"]
        wt.ctx = {"n": n}
        wt.nodes.new(decorated_compare, name="compare1", x="{{n}}", y=orm.Int(limit))
        multiply1 = wt.nodes.new(
            decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
        )
        add1 = wt.nodes.new(decorated_add, name="add1", y=3)
        add1.to_ctx = [["result", "n"]]
        wt.links.new(multiply1.outputs[0], add1.inputs[0])
        return wt

    # -----------------------------------------
    wt = WorkGraph("while")
    add1 = wt.nodes.new(decorated_add, x=orm.Int(25), y=orm.Int(25))
    my_while1 = wt.nodes.new(my_while, n=orm.Int(1))
    add2 = wt.nodes.new(decorated_add, y=orm.Int(2))
    wt.links.new(add1.outputs[0], my_while1.inputs["limit"])
    wt.links.new(my_while1.outputs[0], add2.inputs[0])
    wt.submit(wait=True, timeout=100)
    assert add2.node.outputs.result.value == 63
    assert my_while1.node.outputs.execution_count == 4


def test_while_max_iteration(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    @node.group(outputs=[["ctx.n", "result"]])
    def my_while(n, limit):
        wt = WorkGraph("while_workgraph")
        wt.workgraph_type = "WHILE"
        wt.max_iteration = 3
        wt.conditions = ["compare1.result"]
        wt.ctx = {"n": n}
        wt.nodes.new(decorated_compare, name="compare1", x="{{n}}", y=orm.Int(limit))
        multiply1 = wt.nodes.new(
            decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
        )
        add1 = wt.nodes.new(decorated_add, name="add1", y=3)
        add1.to_ctx = [["result", "n"]]
        wt.links.new(multiply1.outputs[0], add1.inputs[0])
        return wt

    # -----------------------------------------
    wt = WorkGraph("while")
    add1 = wt.nodes.new(decorated_add, x=orm.Int(25), y=orm.Int(25))
    my_while1 = wt.nodes.new(my_while, n=orm.Int(1))
    add2 = wt.nodes.new(decorated_add, y=orm.Int(2))
    wt.links.new(add1.outputs[0], my_while1.inputs["limit"])
    wt.links.new(my_while1.outputs[0], add2.inputs[0])
    wt.submit(wait=True, timeout=100)
    assert add2.node.outputs.result.value < 63
    assert my_while1.node.outputs.execution_count == 3
