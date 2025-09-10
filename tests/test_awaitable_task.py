from aiida_workgraph import WorkGraph, task
import asyncio


def test_awaitable_decorator(decorated_add, capsys):
    @task.awaitable()
    async def awaitable_func(x, y):
        n = 2
        while n > 0:
            n -= 1
            await asyncio.sleep(0.5)
        return x + y

    wg = WorkGraph(name="test_awaitable_decorator")
    awaitable_func1 = wg.add_task(awaitable_func, "awaitable_func1", x=1, y=2)
    add1 = wg.add_task(decorated_add, "add1", x=1, y=awaitable_func1.outputs.result)
    wg.run()
    captured = capsys.readouterr()
    report = captured.out
    assert "Waiting for child processes: awaitable_func1" in report
    assert add1.outputs.result.value == 4
