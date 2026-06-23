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


@task.graph
def add_workflow(x, y) -> int:
    """Async process-type source task (runs as its own sub-process) for use inside a Map zone."""
    return add(x=x, y=y).result


@task()
def maybe_fail(x, y):
    """Add two numbers, but fail for one specific item to exercise the failure path."""
    if x == 1:
        raise ValueError('boom on x==1')
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


def test_map_zone_async_source():
    """Map over an async process-type source task (`@task.graph`) and gather it.

    Regression test for the gather race: when the mapped source is a process-type
    task, the awaitable cascade can reach the gather phase before the gather_item
    clones are scheduled. Before the fix this raised ``KeyError`` and excepted the
    engine; the existing ``test_map_zone`` does not catch it because it maps over
    plain synchronous ``@task`` functions that all succeed.
    """
    x = 1
    n = 3
    with WorkGraph('map_async_source') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = add_workflow(x=map_zone.item.value, y=x).result
            map_zone.gather({'sum1': out1})
        out3 = calc_sum(data=map_zone.outputs.sum1).result
        wg.run()
        # values are 0+1, 1+1, 2+1 -> 1 + 2 + 3 = 6
        assert out3.value == 6


def test_map_zone_failed_source_is_graceful():
    """A single failed mapped task must not except the engine.

    Before the fix, the gather read ``KeyError``-ed on the failed item and the
    engine excepted. Now the gather keeps ``None`` for the failed item and the
    WorkGraph finishes, surfacing the failure via exit status 302.
    """
    n = 3
    with WorkGraph('map_fail') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = maybe_fail(x=map_zone.item.value, y=10).result
            map_zone.gather({'sum1': out1})
        wg.run()
    assert wg.process.exit_status == 302
    assert 'key_1_maybe_fail' in wg.process.exit_message
