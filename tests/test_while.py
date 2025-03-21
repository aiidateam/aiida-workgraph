from aiida_workgraph import WorkGraph, TaskPool


def test_while_instruction(decorated_add, decorated_multiply, decorated_smaller_than):
    from aiida_workgraph import while_

    wg = WorkGraph("test_while")
    wg.ctx = {"n": 1}
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    wg.update_ctx({"n": add1.outputs.result})
    # ---------------------------------------------------------------------
    compare1 = wg.add_task(decorated_smaller_than, name="compare1", x=wg.ctx.n, y=20)
    while_(compare1.outputs["result"], max_iterations=10)(
        wg.add_task(decorated_add, name="add2", x=wg.ctx.n, y=1),
        wg.add_task(
            decorated_multiply,
            name="multiply1",
            x=wg.tasks["add2"].outputs["result"],
            y=2,
        ),
    )
    wg.tasks["add2"].waiting_on.add("add1")
    wg.update_ctx({"n": wg.tasks.multiply1.outputs.result})
    add3 = wg.add_task(decorated_add, name="add3", x=1, y=1)
    wg.add_link(wg.tasks.multiply1.outputs["result"], add3.inputs["x"])
    assert len(wg.tasks) == 6
    assert "while_1" in wg.tasks
    assert len(wg.tasks.while_1.children) == 2
    wg.run()
    assert wg.state == "FINISHED"
    assert wg.tasks.add3.outputs.result.value == 31


def test_while_task(decorated_add, decorated_smaller_than):
    """Test nested while task.
    Also test the max_iteration parameter."""

    def raw_python_code():
        n01 = decorated_add._func(1, 1)
        m = 1
        n = n01
        l1 = 1
        while m < 10:
            n11 = decorated_add._func(1, 1)
            print("n11", n11)
            while n < 5:
                n21 = decorated_add._func(n, n11)
                n22 = decorated_add._func(n21, 1)
                n = n22
                print("n21", n21)
                print("n22", n22)
            niter = 0
            while l1 < 5 and niter < 1:
                n31 = decorated_add._func(l1, 1)
                n32 = decorated_add._func(n31, 1)
                l1 = n32
                niter += 1
                print("n31", n31)
                print("n32", n32)
            n12 = decorated_add._func(m, n32)
            print("n12", n12)
            m = n12

        m = decorated_add._func(m, n31)
        return m

    wg = WorkGraph("test_while_task")
    # set a context variable before running.
    wg.ctx = {
        "m": 1,
        "l": 1,
    }
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    wg.update_ctx({"n": add1.outputs.result})
    # ---------------------------------------------------------------------
    # the `result` of compare1 taskis used as condition
    compare1 = wg.add_task(decorated_smaller_than, name="compare1", x=wg.ctx.m, y=10)
    while1 = wg.add_task(
        TaskPool.workgraph.while_zone, name="while1", conditions=compare1.outputs.result
    )
    add11 = wg.add_task(decorated_add, name="add11", x=1, y=1)
    # ---------------------------------------------------------------------
    compare2 = wg.add_task(decorated_smaller_than, name="compare2", x=wg.ctx.n, y=5)
    while2 = wg.add_task(
        TaskPool.workgraph.while_zone, name="while2", conditions=compare2.outputs.result
    )
    add21 = wg.add_task(decorated_add, name="add21", x=wg.ctx.n, y=add11.outputs.result)
    add21.waiting_on.add("add1")
    add22 = wg.add_task(decorated_add, name="add22", x=add21.outputs.result, y=1)
    wg.update_ctx({"n": add22.outputs.result})
    while2.children.add(["add21", "add22"])
    # ---------------------------------------------------------------------
    compare3 = wg.add_task(decorated_smaller_than, name="compare3", x=wg.ctx.l, y=5)
    while3 = wg.add_task(
        TaskPool.workgraph.while_zone,
        name="while3",
        max_iterations=1,
        conditions=compare3.outputs.result,
    )
    add31 = wg.add_task(decorated_add, name="add31", x=wg.ctx.l, y=1)
    add31.waiting_on.add("add22")
    add32 = wg.add_task(decorated_add, name="add32", x=add31.outputs.result, y=1)
    wg.update_ctx({"l": add32.outputs.result})
    while3.children.add(["add31", "add32"])
    # ---------------------------------------------------------------------
    add12 = wg.add_task(decorated_add, name="add12", x=wg.ctx.m, y=add32.outputs.result)
    wg.update_ctx({"m": add12.outputs.result})
    while1.children.add(["add11", "while2", "while3", "add12", "compare2", "compare3"])
    # ---------------------------------------------------------------------
    add2 = wg.add_task(
        decorated_add, name="add2", x=add12.outputs.result, y=add31.outputs.result
    )
    # wg.submit(wait=True, timeout=100)
    wg.run()
    # print out the node labels and the results for debugging
    # print("link node label, result")
    # for link in wg.process.base.links.get_outgoing().all():
    #     if isinstance(link.node, orm.ProcessNode):
    #         print(link.node.label, link.node.outputs.result)
    # print("assert check")
    # assert add2.outputs.result.value.value == raw_python_code().value
    assert add2.outputs.result.value.value == 18
