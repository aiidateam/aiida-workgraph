from aiida_workgraph import (
    WorkGraph,
    task,
    active_graph,
    active_map_zone,
    active_if_zone,
)
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


def test_map_instruction(add_code):
    x = 1
    y = 2
    n = 3
    with active_graph(WorkGraph("add_graph")) as wg:
        wg.add_task(generate_data, name="generate_data", n=n)
        with active_map_zone(wg.tasks.generate_data.outputs.result) as map_zone:
            map_zone.add_task(
                ArithmeticAddCalculation,
                name="add1",
                x=map_zone.item,
                y=y,
                code=add_code,
            )
            with active_if_zone(wg.tasks.add1.outputs.sum < 4) as if_zone1:
                if_zone1.add_task(
                    ArithmeticAddCalculation,
                    name="add2",
                    x=x,
                    y=wg.tasks.add1.outputs.sum,
                    code=add_code,
                )
                wg.update_ctx({"sum": wg.tasks.add2.outputs.sum})
            with active_if_zone(wg.tasks.add1.outputs.sum >= 4) as if_zone2:
                if_zone2.add_task(
                    ArithmeticAddCalculation,
                    name="add3",
                    x=3,
                    y=wg.tasks.add1.outputs.sum,
                    code=add_code,
                )
                wg.update_ctx({"sum": wg.tasks.add3.outputs.sum})
        wg.add_task(calc_sum, name="calc_sum1", kwargs=wg.ctx.sum)
        wg.tasks.calc_sum1.waiting_on.add([wg.tasks.add2, wg.tasks.add3])
        wg.run()
        assert wg.tasks.calc_sum1.outputs.sum.value == 7
