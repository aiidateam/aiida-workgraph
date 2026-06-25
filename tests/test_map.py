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


@task()
def echo(text: str) -> str:
    """Return the input unchanged."""
    return text


@task()
def join_keys(data: Annotated[dict, dynamic(orm.Str)]) -> str:
    """Join all provided strings in sorted order."""
    return ','.join(sorted(data.values()))


def test_map_zone():
    x = 1
    y = 2
    n = 3
    with WorkGraph('add_graph') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = add(x=map_zone.value, y=x).result
            out2 = add(x=map_zone.value, y=y).result
            map_zone.gather({'sum1': out1, 'sum2': out2})
        out3 = calc_sum(data=map_zone.outputs.sum1).result
        out4 = calc_sum(data=map_zone.outputs.sum2).result
        wg.run()
        assert out3.value == 6
        assert out4.value == 9


def test_map_value_and_key():
    """z.value and z.key resolve to the same per-element item and both flow (#785)."""
    with WorkGraph('map_value_key') as wg:
        data = generate_data(n=2).data
        with Map(data) as z:
            s = add(x=z.value, y=10).result
            k = echo(text=z.key).result
            z.gather({'sum': s, 'k': k})
        total = calc_sum(data=z.outputs.sum).result
        joined = join_keys(data=z.outputs.k).result
        wg.run()
    # values 0, 1 -> (0+10)+(1+10) = 21; keys are the user's own source keys
    assert total.value == 21
    assert joined.value == 'key_0,key_1'
