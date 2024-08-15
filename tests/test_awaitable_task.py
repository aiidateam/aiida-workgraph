from aiida_workgraph import WorkGraph, task
import asyncio
from aiida.cmdline.utils.common import get_workchain_report


def test_awaitable_task(decorated_add):
    @task.awaitable()
    async def awaitable_func(x, y):
        n = 2
        while n > 0:
            n -= 1
            await asyncio.sleep(0.5)
        return x + y

    wg = WorkGraph(name="test_awaitable")
    awaitable1 = wg.add_task(awaitable_func, "awaitable_func1", x=1, y=2)
    add1 = wg.add_task(decorated_add, "add1", x=1, y=awaitable1.outputs["result"])
    wg.run()
    report = get_workchain_report(wg.process, "REPORT")
    assert "Waiting for child processes: awaitable_func1" in report
    assert add1.outputs["result"].value == 4
