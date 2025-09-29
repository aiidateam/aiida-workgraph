from aiida_workgraph import WorkGraph, task
import asyncio
from aiida.cmdline.utils.common import get_workchain_report


def test_async_decorator(decorated_add):
    @task()
    async def async_func(x, y):
        n = 2
        while n > 0:
            n -= 1
            await asyncio.sleep(0.5)
        return x + y

    wg = WorkGraph(name='test_async_functioin')
    async_func1 = wg.add_task(async_func, 'async_func1', x=1, y=2)
    add1 = wg.add_task(decorated_add, 'add1', x=1, y=async_func1.outputs.result)
    wg.run()
    assert add1.outputs.result.value == 4
    report = get_workchain_report(wg.process, 'REPORT')
    assert 'Waiting for child processes: ' in report
