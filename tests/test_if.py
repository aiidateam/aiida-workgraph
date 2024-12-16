from aiida_workgraph import WorkGraph


def test_if_task(decorated_add, decorated_multiply, decorated_compare):
    """Test the If task."""

    wg = WorkGraph("test_if")
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    condition1 = wg.add_task(decorated_compare, name="condition1", x=1, y=0)
    if_zone = wg.add_task("If", name="if_true", conditions=condition1.outputs.result)
    add2 = if_zone.add_task(decorated_add, name="add2", x=add1.outputs.result, y=2)
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x=add1.outputs.result, y=2
    )
    if2 = wg.add_task(
        "If",
        name="if_false",
        conditions=condition1.outputs.result,
        invert_condition=True,
    )
    if2.children.add("multiply1")
    # ---------------------------------------------------------------------
    select1 = wg.add_task(
        "workgraph.select",
        name="select1",
        true=add2.outputs.result,
        false=multiply1.outputs.result,
        condition=condition1.outputs.result,
    )
    add3 = wg.add_task(decorated_add, name="add3", x=select1.outputs.result, y=1)
    wg.run()
    assert add3.outputs.result.value == 5


def test_empty_if_task():
    """Test the If task with no children."""

    wg = WorkGraph("test_empty_if")
    wg.add_task("If", name="if_true")
    wg.run()
    assert wg.state == "FINISHED"
