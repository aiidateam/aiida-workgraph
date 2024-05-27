from aiida_workgraph import node, WorkGraph
from aiida import load_profile, orm

load_profile()


def test_while(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    wg = WorkGraph("while_workgraph")
    wg.workgraph_type = "WHILE"
    wg.conditions = ["compare1.result"]
    wg.context = {"n": 1}
    wg.max_iteration = 10
    wg.nodes.new(decorated_compare, name="compare1", x="{{n}}", y=50)
    multiply1 = wg.nodes.new(
        decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
    )
    add1 = wg.nodes.new(decorated_add, name="add1", y=3)
    add1.to_context = [["result", "n"]]
    wg.links.new(multiply1.outputs["result"], add1.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert wg.execution_count == 4
    assert wg.nodes["add1"].outputs["result"].value == 61


def test_while_graph_builder(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    @node.graph_builder(outputs=[["context.n", "result"]])
    def my_while(n=0, limit=100):
        wg = WorkGraph("while_workgraph")
        wg.workgraph_type = "WHILE"
        wg.conditions = ["compare1.result"]
        wg.context = {"n": n}
        wg.max_iteration = 10
        wg.nodes.new(decorated_compare, name="compare1", x="{{n}}", y=orm.Int(limit))
        multiply1 = wg.nodes.new(
            decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
        )
        add1 = wg.nodes.new(decorated_add, name="add1", y=3)
        add1.to_context = [["result", "n"]]
        wg.links.new(multiply1.outputs["result"], add1.inputs["x"])
        return wg

    # -----------------------------------------
    wg = WorkGraph("while")
    add1 = wg.nodes.new(decorated_add, name="add1", x=orm.Int(25), y=orm.Int(25))
    my_while1 = wg.nodes.new(my_while, n=orm.Int(1))
    add2 = wg.nodes.new(decorated_add, name="add2", y=orm.Int(2))
    wg.links.new(add1.outputs["result"], my_while1.inputs["limit"])
    wg.links.new(my_while1.outputs["result"], add2.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert add2.outputs["result"].value == 63
    assert my_while1.node.outputs.execution_count == 4
    assert my_while1.outputs["result"].value == 61


def test_while_max_iteration(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    @node.graph_builder(outputs=[["context.n", "result"]])
    def my_while(n=0, limit=100):
        wg = WorkGraph("while_workgraph")
        wg.workgraph_type = "WHILE"
        wg.max_iteration = 3
        wg.conditions = ["compare1.result"]
        wg.context = {"n": n}
        wg.nodes.new(decorated_compare, name="compare1", x="{{n}}", y=orm.Int(limit))
        multiply1 = wg.nodes.new(
            decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
        )
        add1 = wg.nodes.new(decorated_add, name="add1", y=3)
        add1.to_context = [["result", "n"]]
        wg.links.new(multiply1.outputs["result"], add1.inputs["x"])
        return wg

    # -----------------------------------------
    wg = WorkGraph("while")
    add1 = wg.nodes.new(decorated_add, name="add1", x=orm.Int(25), y=orm.Int(25))
    my_while1 = wg.nodes.new(my_while, n=orm.Int(1))
    add2 = wg.nodes.new(decorated_add, name="add2", y=orm.Int(2))
    wg.links.new(add1.outputs["result"], my_while1.inputs["limit"])
    wg.links.new(my_while1.outputs["result"], add2.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert add2.outputs["result"].value < 63
    assert my_while1.node.outputs.execution_count == 3
