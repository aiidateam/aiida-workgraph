from aiida_workgraph import WorkGraph, task, map_, if_
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida import orm


@task.calcfunction(outputs=["result"])
def generate_data(n: int) -> dict:
    """Generate a dictionary of integers."""
    result = {f"key_{i}": orm.Int(i) for i in range(n.value)}
    return {"result": result}


@task.calcfunction(outputs=["sum"])
def calc_sum(**kwargs) -> dict:
    """Compute the sum of all provided values."""
    return {"sum": sum(kwargs.values())}


def test_map_instruction(add_code, decorated_smaller_than):
    x = 1
    y = 2
    n = 3
    wg = WorkGraph("add_graph")
    wg.add_task(generate_data, name="generate_data", n=n)
    map_(wg.tasks.generate_data.outputs.result)(
        wg.add_task(
            ArithmeticAddCalculation,
            name="add1",
            x=map_.default_placeholder,
            y=y,
            code=add_code,
        ),
        wg.add_task(
            decorated_smaller_than,
            name="smaller_than",
            x=wg.tasks.add1.outputs.sum,
            y=4,
        ),
        if_(wg.tasks.smaller_than.outputs.result)(
            wg.add_task(
                ArithmeticAddCalculation,
                name="add2",
                x=x,
                y=wg.tasks.add1.outputs.sum,
                code=add_code,
            ),
            wg.tasks.add2.set_context({"sum": "sum"}),
        ).else_(
            wg.add_task(
                ArithmeticAddCalculation,
                name="add3",
                x=3,
                y=wg.tasks.add1.outputs.sum,
                code=add_code,
            ),
            wg.tasks.add2.set_context({"sum": "sum"}),
        ),
    )
    wg.add_task(calc_sum, name="calc_sum1", kwargs=wg.tasks.add2.outputs.sum)
    wg.run()
    assert wg.tasks.calc_sum1.outputs.sum.value == 7
