from aiida_workgraph import WorkGraph, task


def test_create_task_from_calcJob(add_code) -> None:
    """Test creating a task from a CalcJob."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    AddTask = task()(ArithmeticAddCalculation)
    with WorkGraph() as wg:
        outputs1 = AddTask(x=2, y=3, code=add_code)
        outputs2 = AddTask(x=outputs1.sum, y=3, code=add_code)
        wg.run()
    assert outputs2.sum.value == 8
    assert outputs1._node._spec.mode == "decorator_build"
    assert outputs1._node.get_executor().callable == ArithmeticAddCalculation
