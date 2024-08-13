import pytest
from aiida_workgraph import task, WorkGraph
from aiida import orm


@pytest.mark.usefixtures("started_daemon_client")
def test_while_task(decorated_add, decorated_compare):
    """Test nested while task.
    Also test the max_iteration parameter."""

    def raw_python_code():
        n01 = decorated_add(1, 1)
        m = 1
        n = n01
        l1 = 1
        while m < 10:
            n11 = decorated_add(1, 1)
            print("n11", n11)
            while n < 5:
                n21 = decorated_add(n, n11)
                n22 = decorated_add(n21, 1)
                n = n22
                print("n21", n21)
                print("n22", n22)
            niter = 0
            while l1 < 5 and niter < 1:
                n31 = decorated_add(l1, 1)
                n32 = decorated_add(n31, 1)
                l1 = n32
                niter += 1
                print("n31", n31)
                print("n32", n32)
            n12 = decorated_add(m, n32)
            print("n12", n12)
            m = n12

        m = decorated_add(m, n31)
        return m

    wg = WorkGraph("test_while_task")
    # set a context variable before running.
    wg.context = {
        "should_run1": True,
        "should_run2": True,
        "should_run3": True,
        "m": 1,
        "l": 1,
    }
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    add1.set_context({"result": "n"})
    # ---------------------------------------------------------------------
    while1 = wg.add_task("While", name="while1", conditions=["should_run1"])
    add11 = wg.add_task(decorated_add, name="add11", x=1, y=1)
    # ---------------------------------------------------------------------
    while2 = wg.add_task("While", name="while2", conditions=["should_run2"])
    add21 = wg.add_task(
        decorated_add, name="add21", x="{{n}}", y=add11.outputs["result"]
    )
    add21.wait.append("add1")
    add22 = wg.add_task(decorated_add, name="add22", x=add21.outputs["result"], y=1)
    add22.set_context({"result": "n"})
    compare2 = wg.add_task(
        decorated_compare, name="compare2", x=add22.outputs["result"], y=5
    )
    compare2.set_context({"result": "should_run2"})
    while2.children = ["add21", "add22", "compare2"]
    # ---------------------------------------------------------------------
    while3 = wg.add_task(
        "While", name="while3", max_iterations=1, conditions=["should_run3"]
    )
    add31 = wg.add_task(decorated_add, name="add31", x="{{l}}", y=1)
    add31.wait.append("add22")
    add32 = wg.add_task(decorated_add, name="add32", x=add31.outputs["result"], y=1)
    add32.set_context({"result": "l"})
    compare3 = wg.add_task(
        decorated_compare, name="compare3", x=add32.outputs["result"], y=5
    )
    compare3.set_context({"result": "should_run3"})
    while3.children = ["add31", "add32", "compare3"]
    # ---------------------------------------------------------------------
    add12 = wg.add_task(
        decorated_add, name="add12", x="{{m}}", y=add32.outputs["result"]
    )
    add12.set_context({"result": "m"})
    compare1 = wg.add_task(
        decorated_compare, name="compare1", x=add12.outputs["result"], y=10
    )
    compare1.set_context({"result": "should_run1"})
    while1.children = ["add11", "while2", "while3", "add12", "compare1"]
    # the `result` of compare1 taskis used as condition
    # ---------------------------------------------------------------------
    add2 = wg.add_task(
        decorated_add, name="add2", x=add12.outputs["result"], y=add31.outputs["result"]
    )
    # wg.submit(wait=True, timeout=100)
    wg.run()
    assert add2.outputs["result"].value == raw_python_code()


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
