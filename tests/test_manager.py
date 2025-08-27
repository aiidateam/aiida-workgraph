from aiida_workgraph import WorkGraph, task
from aiida import orm
from aiida_workgraph.manager import (
    If,
    While,
    Map,
)
from aiida_workgraph import socket_spec as spec


@task()
def is_even(x: orm.Int):
    return x % 2 == 0


def test_while_and_if(decorated_add):
    """Sum of even numbers up to N."""

    N = 5
    with WorkGraph() as wg:
        wg.ctx = {"n": 1, "total": 0}
        with While(wg.ctx.n < N):
            with If(wg.ctx.n % 2 == 0):
                outputs = decorated_add(x=wg.ctx.total, y=wg.ctx.n)
                wg.ctx.total = outputs.result
            n = wg.ctx.n + 1
            wg.ctx.n = n
            n << outputs.result
        wg.run()

    assert outputs.result.value == 6


def test_map(decorated_add):
    """"""

    @task
    def generate_list(N) -> spec.namespace(result=spec.dynamic(any)):
        """Generate a list of N items."""
        return {"result": {f"item_{i}": i for i in range(1, N + 1)}}

    @task()
    def sum_values(**items):
        return sum(items.values())

    N = 5
    with WorkGraph() as wg:
        outputs = generate_list(N)
        with Map(source_socket=outputs.result) as map_zone:
            outputs1 = decorated_add(x=map_zone.item, y=1)
        outputs2 = sum_values(items=outputs1.result)
        wg.run()
        # wg.to_html("test_map.html")
    assert outputs2.result.value == 20
