from aiida_workgraph import task, WorkGraph


def test_if_task(decorated_add, decorated_multiply, decorated_compare):
    """Test the If task."""

    @task()
    def select(condition, x, y):
        """Output the result based on the condition."""
        if condition:
            return x
        return y

    wg = WorkGraph("test_if")
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    condition1 = wg.add_task(decorated_compare, name="condition1", x=1, y=0)
    add2 = wg.add_task(decorated_add, name="add2", x=add1.outputs["result"], y=2)
    if1 = wg.add_task("If", name="if_true", conditions=condition1.outputs["result"])
    if1.children.append("add2")
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x=add1.outputs["result"], y=2
    )
    if2 = wg.add_task(
        "If",
        name="if_false",
        conditions=condition1.outputs["result"],
        invert_condition=True,
    )
    if2.children.append("multiply1")
    # ---------------------------------------------------------------------
    select1 = wg.add_task(
        select,
        name="select1",
        x=add2.outputs["result"],
        y=multiply1.outputs["result"],
        condition=condition1.outputs["result"],
    )
    add3 = wg.add_task(decorated_add, name="add3", x=select1.outputs["result"], y=1)
    wg.run()
    assert add3.outputs["result"].value == 5
