from aiida_workgraph import WorkGraph, task
from aiida import orm
from aiida_workgraph.manager import (
    active_graph,
    active_if_zone,
    active_while_zone,
    active_map_zone,
)


@task()
def is_even(x: orm.Int):
    return x % 2 == 0


def test_while_and_if(decorated_add):
    """Sum of even numbers up to N."""

    N = 5
    with active_graph(WorkGraph()) as wg:
        wg.ctx = {"n": 1, "total": 0}
        with active_while_zone(wg.ctx.n < N):
            with active_if_zone(wg.ctx.n % 2 == 0):
                total = decorated_add(x=wg.ctx.total, y=wg.ctx.n)
                wg.update_ctx({"total": total})
            n = wg.ctx.n + 1
            wg.update_ctx({"n": n})
            n._waiting_on.add(total)
        wg.run()

    assert total.value == 6


def test_map(decorated_add):
    """"""

    @task()
    def generate_list(N):
        return {"result": {f"item_{i}": i for i in range(1, N + 1)}}

    @task()
    def sum_values(**items):
        return sum(items.values())

    N = 5
    with active_graph(WorkGraph()) as wg:
        result = generate_list(N)
        with active_map_zone(source_socket=result) as map_zone:
            result = decorated_add(x=map_zone.item, y=1)
        total = sum_values(items=result)
        total._waiting_on.add(result)
        wg.run()
        # wg.to_html("test_map.html")
    assert total.value == 20
