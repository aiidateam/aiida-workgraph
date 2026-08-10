import pytest
from aiida_workgraph import (
    WorkGraph,
    task,
    Map,
    If,
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
def is_small(x) -> bool:
    """True for x < 2, to drive an If branch per map item."""
    return x < 2


@task()
def double(x):
    """Double a number."""
    return x * 2


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
            out1 = add_workflow(x=map_zone.value, y=x).result
            map_zone.gather({'sum1': out1})
        out3 = calc_sum(data=map_zone.outputs.sum1).result
        wg.run()
        # values are 0+1, 1+1, 2+1 -> 1 + 2 + 3 = 6
        assert out3.value == 6


def test_map_zone_failed_iteration_fails_the_zone():
    """A failed mapped iteration must fail the zone, not except and not under-report.

    Before the fix the gather ``KeyError``-ed on the item that produced no
    result, excepting the engine and burying the real cause. The zone now goes
    FAILED and gathers nothing, rather than reporting FINISHED with a namespace
    that is silently missing the failed item.
    """
    n = 3
    with WorkGraph('map_fail') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = maybe_fail(x=map_zone.value, y=10).result
            map_zone.gather({'sum1': out1})
        wg.run()
    assert map_zone.state == 'FAILED'
    assert wg.process.exit_status == 302
    assert 'key_1_maybe_fail' in wg.process.exit_message


def test_map_zone_failed_iteration_skips_downstream():
    """A failed zone must SKIP its downstream tasks, like any ordinary failure.

    Otherwise the consumer runs on the missing gather output, fails on its own,
    and pollutes the report with a failure that is only a consequence of the
    zone's.
    """
    n = 3
    with WorkGraph('map_fail_downstream') as wg:
        data = generate_data(n=n).data
        with Map(data) as map_zone:
            out1 = maybe_fail(x=map_zone.value, y=10).result
            map_zone.gather({'sum1': out1})
        calc_sum(data=map_zone.outputs.sum1)
        wg.run()
    assert map_zone.state == 'FAILED'
    assert wg.process.get_task_state('calc_sum') == 'SKIPPED'
    assert 'calc_sum' not in wg.process.exit_message


def test_map_gather_rejects_outside_zone_source():
    """A gather source produced outside the zone fails loudly at build time.

    It has no per-item clones, so it would gather to an empty namespace; reject
    it with a clear message instead of silently under-reporting.
    """
    with pytest.raises(ValueError, match='outside the Map zone'):
        with WorkGraph('map_reject_outside'):
            data = generate_data(n=2).data
            outside = add(x=1, y=2).result  # produced OUTSIDE the Map zone
            with Map(data) as map_zone:
                inner = add(x=map_zone.value, y=0).result
                map_zone.gather({'inner': inner, 'outside': outside})


def test_map_if_branch_does_not_fail_zone():
    """A deliberately-untaken If branch inside a Map is a skip, not a failure.

    For items where the condition is false the branch (and gather source) is
    SKIPPED; that item gathers None and the zone still finishes, rather than the
    zone failing as it would for a genuine error.
    """
    with WorkGraph('map_if') as wg:
        data = generate_data(n=3).data  # values 0, 1, 2
        with Map(data) as map_zone:
            with If(is_small(x=map_zone.value).result):  # False for value 2
                out = double(x=map_zone.value).result
            map_zone.gather({'doubled': out})
        wg.run()
    assert map_zone.state == 'FINISHED'
    assert wg.process.exit_status == 0
    assert wg.process.get_task_state('key_2_double') == 'SKIPPED'


def test_map_upstream_failure_fails_zone():
    """An upstream error that only SKIPs the gather source still fails the zone.

    Distinguishes a real failure from a deliberate skip: the gather source here
    is SKIPPED (its input errored), but a FAILED clone exists in that iteration,
    so the zone fails rather than gathering None and finishing.
    """
    with WorkGraph('map_upstream_fail') as wg:
        data = generate_data(n=3).data
        with Map(data) as map_zone:
            failed = maybe_fail(x=map_zone.value, y=0).result  # FAILS for key_1
            out = add(x=failed, y=100).result  # gather source; SKIPPED for key_1
            map_zone.gather({'out': out})
        wg.run()
    assert map_zone.state == 'FAILED'
    assert wg.process.get_task_state('key_1_add') == 'SKIPPED'
    assert 'key_1' in wg.process.exit_message


def test_map_unrelated_body_failure_does_not_fail_zone():
    """An un-gathered body task failing must not fail the zone or discard the gather.

    Only the gather sources decide completeness; an unrelated failure surfaces via the
    ordinary exit-302 path, and the (complete) gathered output and its consumer are
    unaffected.
    """
    with WorkGraph('map_unrelated_fail') as wg:
        data = generate_data(n=3).data
        with Map(data) as map_zone:
            good = add(x=map_zone.value, y=1).result
            maybe_fail(x=map_zone.value, y=0)  # FAILS for key_1; its output is not gathered
            map_zone.gather({'a': good})
        total = calc_sum(data=map_zone.outputs.a).result
        wg.run()
    assert map_zone.state == 'FINISHED'
    assert wg.process.get_task_state('calc_sum') == 'FINISHED'
    assert wg.process.exit_status == 302  # only the un-gathered maybe_fail failed
    assert 'map_zone' not in wg.process.exit_message
    assert total.value == 6  # gather complete: (0+1)+(1+1)+(2+1)


def test_map_info_keeps_gather_edges():
    """`gather_item` is not cloned, but its source edges stay in `map_info` for the GUI.

    They use template names (like every other map_info link) and resolve to the listed
    children, so the web UI can still draw the per-item gather.
    """
    with WorkGraph('map_info') as wg:
        data = generate_data(n=2).data
        with Map(data) as map_zone:
            out = add(x=map_zone.value, y=1).result
            map_zone.gather({'a': out})
        wg.run()
    mi = wg.process.get_task_map_info('map_zone')
    assert 'gather_item' in mi['children']
    gather_edges = [link for link in mi['links'] if link['to_task'] == 'gather_item']
    assert [e['from_task'] for e in gather_edges] == ['add']
    # no dangling references: every link's endpoints are listed children
    assert all(link['from_task'] in mi['children'] and link['to_task'] in mi['children'] for link in mi['links'])


def test_map_gather_rejects_atomically():
    """A rejected gather() adds no specs, so the zone stays rebuildable.

    Otherwise a caught error would leave a half-built zone whose sockets already
    exist on the retry.
    """
    with WorkGraph('map_reject_atomic') as wg:
        data = generate_data(n=2).data
        outside = add(x=1, y=2).result  # outside the zone
        with Map(data) as map_zone:
            good = add(x=map_zone.value, y=0).result
            with pytest.raises(ValueError, match='outside the Map zone'):
                map_zone.gather({'a': good, 'ext': outside})
            # the rejected call left nothing behind, so a valid gather still works
            map_zone.gather({'a': good})
        out = calc_sum(data=map_zone.outputs.a).result
        wg.run()
    assert out.value == 1  # (0+0) + (1+0)
