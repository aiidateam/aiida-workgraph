import pytest
from aiida_workgraph import task, WorkGraph
from aiida import orm


@pytest.mark.usefixtures("started_daemon_client")
def test_while_task(decorated_add, decorated_multiply, decorated_compare):
    wg = WorkGraph("test_while_task")
    # set a context variable before running.
    wg.context = {"should_run": True}
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    add1.set_context({"result": "n"})
    # ---------------------------------------------------------------------
    add2 = wg.add_task(decorated_add, name="add2", x="{{n}}", y=1)
    add2.wait.append("add1")
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x=add2.outputs["result"], y=2
    )
    # update the context variable
    multiply1.set_context({"result": "n"})
    compare1 = wg.add_task(
        decorated_compare, name="compare1", x=multiply1.outputs["result"], y=30
    )
    compare1.set_context({"result": "should_run"})
    wg.add_task(
        "While",
        max_iterations=100,
        conditions=["should_run"],
        tasks=["add2", "multiply1", "compare1"],
    )
    # the `result` of compare1 taskis used as condition
    # ---------------------------------------------------------------------
    add3 = wg.add_task(decorated_add, name="add3", x=1, y=1)
    wg.add_link(multiply1.outputs["result"], add3.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert wg.tasks["add3"].outputs["result"].value == 31


def test_while_workgraph(decorated_add, decorated_multiply, decorated_compare):
    # Create a WorkGraph will repeat itself based on the conditions
    wg = WorkGraph("while_workgraph")
    wg.workgraph_type = "WHILE"
    wg.conditions = ["compare1.result"]
    wg.context = {"n": 1}
    wg.max_iteration = 10
    wg.add_task(decorated_compare, name="compare1", x="{{n}}", y=20)
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
    )
    add1 = wg.add_task(decorated_add, name="add1", y=3)
    add1.set_context({"result": "n"})
    wg.add_link(multiply1.outputs["result"], add1.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert wg.execution_count == 3
    assert wg.tasks["add1"].outputs["result"].value == 29


@pytest.mark.usefixtures("started_daemon_client")
def test_while_graph_builder(decorated_add, decorated_multiply, decorated_compare):
    """Test the while WorkGraph in graph builder.
    Also test the max_iteration parameter."""

    @task.graph_builder(outputs=[{"name": "result", "from": "context.n"}])
    def my_while(n=0, limit=100):
        wg = WorkGraph("while_workgraph")
        wg.workgraph_type = "WHILE"
        wg.max_iteration = 2
        wg.conditions = ["compare1.result"]
        wg.context = {"n": n}
        wg.add_task(decorated_compare, name="compare1", x="{{n}}", y=orm.Int(limit))
        multiply1 = wg.add_task(
            decorated_multiply, name="multiply1", x="{{ n }}", y=orm.Int(2)
        )
        add1 = wg.add_task(decorated_add, name="add1", y=3)
        add1.set_context({"result": "n"})
        wg.add_link(multiply1.outputs["result"], add1.inputs["x"])
        return wg

    # -----------------------------------------
    wg = WorkGraph("while")
    add1 = wg.add_task(decorated_add, name="add1", x=orm.Int(25), y=orm.Int(25))
    my_while1 = wg.add_task(my_while, n=orm.Int(1))
    add2 = wg.add_task(decorated_add, name="add2", y=orm.Int(2))
    wg.add_link(add1.outputs["result"], my_while1.inputs["limit"])
    wg.add_link(my_while1.outputs["result"], add2.inputs["x"])
    wg.submit(wait=True, timeout=100)
    assert add2.outputs["result"].value < 31
    assert my_while1.node.outputs.execution_count == 2
