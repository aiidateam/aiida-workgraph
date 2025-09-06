from aiida_workgraph import WorkGraph, if_, TaskPool


def test_if_instruction(decorated_add, decorated_multiply, decorated_smaller_than):
    """Test the If instruction."""

    wg = WorkGraph()
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=-2)
    condition1 = wg.add_task(
        decorated_smaller_than, name="condition1", x=add1.outputs.result, y=0
    )
    if_(condition1.outputs.result)(
        wg.add_task(decorated_add, name="add2", x=add1.outputs.result, y=2),
    ).else_(
        wg.add_task(decorated_multiply, name="multiply1", x=add1.outputs.result, y=2),
    )
    assert len(wg.tasks) == 9
    assert "if_true_1" in wg.tasks
    assert "if_false_2" in wg.tasks
    assert len(wg.tasks["if_true_1"].children) == 1
    assert len(wg.tasks["if_false_2"].children) == 1
    wg.run()
    assert wg.state == "FINISHED"
    assert wg.tasks["add2"].outputs.result.value == 1
    assert wg.tasks["multiply1"].outputs.result.value is None


def test_if_task(decorated_add, decorated_multiply, decorated_smaller_than):
    """Test the If task."""

    wg = WorkGraph("test_if")
    add1 = wg.add_task(decorated_add, name="add1", x=1, y=1)
    condition1 = wg.add_task(decorated_smaller_than, name="condition1", x=1, y=0)
    if_zone = wg.add_task(
        TaskPool.workgraph.if_zone, name="if_true", conditions=condition1.outputs.result
    )
    add2 = if_zone.add_task(decorated_add, name="add2", x=add1.outputs.result, y=2)
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x=add1.outputs.result, y=2
    )
    if2 = wg.add_task(
        TaskPool.workgraph.if_zone,
        name="if_false",
        conditions=condition1.outputs.result,
        invert_condition=True,
    )
    if2.children.add("multiply1")
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
    wg.add_task(TaskPool.workgraph.if_zone, name="if_true")
    wg.run()
    assert wg.state == "FINISHED"
