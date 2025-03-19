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


def test_while_and_if(decorated_smaller_than, decorated_add):
    """Sum of even numbers up to N."""

    N = 5
    with active_graph(WorkGraph()) as wg:
        wg.ctx = {"n": 1, "total": 0}
        result = decorated_smaller_than(x=wg.ctx.n, y=N)
        with active_while_zone(result):
            result = is_even(x=wg.ctx.n)
            with active_if_zone(result):
                total = decorated_add(x=wg.ctx.total, y=wg.ctx.n)
                wg.update_ctx({"total": total})
            n = decorated_add(x=wg.ctx.n, y=1)
            wg.update_ctx({"n": n})
            n._node.waiting_on.add(total._node)
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
        with active_map_zone(source_socket=result, placeholder="item"):
            result = decorated_add(x="{{item}}", y=1)
            print("result: ", result)
        total = sum_values(items=result)
        total._node.waiting_on.add(result._node)
        wg.run()
        # wg.to_html("test_map.html")
    assert total.value == 20
