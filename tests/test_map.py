from aiida_workgraph import (
    WorkGraph,
    task,
    Map,
    namespace,
    dynamic,
)
from aiida import orm
from typing import Annotated


@task()
def generate_data(n: int) -> Annotated[dict, namespace(data=dynamic(int))]:
    """Generate a dictionary of integers."""
    result = {f'key_{i}': i for i in range(n)}
    return {'data': result}


@task()
def add(x, y):
    """Add two numbers."""
    return x + y


@task()
def calc_sum(data: Annotated[dict, dynamic(orm.Int)]) -> float:
    """Compute the sum of all provided values."""
    return sum(data.values())


def test_map_zone():
    x = 1
    y = 2
    n = 3
    with WorkGraph('add_graph') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = add(x=map_zone.item.value, y=x).result
            out2 = add(x=map_zone.item.value, y=y).result
            map_zone.gather({'sum1': out1, 'sum2': out2})
        out3 = calc_sum(data=map_zone.outputs.sum1).result
        out4 = calc_sum(data=map_zone.outputs.sum2).result
        wg.run()
        assert out3.value == 6
        assert out4.value == 9
