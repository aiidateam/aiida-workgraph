import pytest
from aiida_workgraph import WorkGraph
from typing import Callable
from aiida.engine import run, submit


def test_normal_function_run(decorated_normal_add: Callable, decorated_add: Callable) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name='test_normal_function_run')
    add1 = wg.add_task(decorated_normal_add, 'add1', x=2, y=3)
    add2 = wg.add_task(decorated_add, 'add2', x=6)
    wg.add_link(add1.outputs.result, add2.inputs['y'])
    run(wg)
    assert wg.tasks.add2.outputs.result.value == 11


@pytest.mark.usefixtures('started_daemon_client')
def test_normal_function_submit(decorated_normal_add: Callable, decorated_add: Callable) -> None:
    """Run simple calcfunction."""
    wg = WorkGraph(name='test_normal_function_submit')
    add1 = wg.add_task(decorated_normal_add, 'add1', x=2, y=3)
    add2 = wg.add_task(decorated_add, 'add2', x=6)
    wg.add_link(add1.outputs.result, add2.inputs['y'])
    submit(wg)
    wg.wait(timeout=120)
    assert wg.tasks.add2.outputs.result.value == 11
