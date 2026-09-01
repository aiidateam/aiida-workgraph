import pytest
from typing import Callable


@pytest.mark.usefixtures('started_daemon_client')
def test_failed_node(decorated_sqrt: Callable, decorated_add: Callable) -> None:
    """Submit simple calcfunction."""
    from aiida_workgraph import WorkGraph
    from aiida.orm import Float

    wg = WorkGraph(name='test_failed_node')
    wg.add_task(decorated_add, 'add1', x=Float(1), y=Float(2))
    sqrt1 = wg.add_task(decorated_sqrt, 'sqrt1', x=Float(-1))
    wg.add_task(decorated_sqrt, 'sqrt2', x=sqrt1.outputs.result)
    wg.submit(wait=True)
    # print("results: ", results[])
    assert wg.process.exit_status == 302
    assert (
        wg.process.exit_message
        == "WorkGraph finished, but tasks: ['sqrt1'] failed. Thus all their child tasks are skipped."
    )


@pytest.mark.parametrize('decorator_name', ['calcfunction', 'workfunction'])
def test_get_task_process_for_excepted_function_task(decorator_name: str) -> None:
    """A raising calcfunction/workfunction still leaves its process reference retrievable.

    Regression test for aiidateam/aiida-workgraph#809: ``run_get_node`` stores the
    excepted node and CALL-links it before re-raising the wrapped function's
    exception, but the engine used to lose that reference on the way out, so
    ``get_task_process`` returned ``None`` for a task that provenance shows did run.
    """
    from aiida.orm import Int, ProcessNode

    from aiida_workgraph import WorkGraph, task
    from aiida_workgraph.orm.utils import deserialize_safe

    decorator = getattr(task, decorator_name)

    @decorator()
    def fails(x):
        raise ValueError('deliberate failure')
        return x  # unreachable; gives the task's outputs a 'result' socket to link from

    @decorator()
    def downstream(x):
        return x

    wg = WorkGraph(name=f'test_excepted_{decorator_name}_process_reference')
    fail_task = wg.add_task(fails, 'fail1', x=Int(1))
    wg.add_task(downstream, 'after1', x=fail_task.outputs.result)
    wg.run()

    node = wg.process
    assert node.exit_status == 302
    assert node.is_finished_ok is False
    assert node.get_task_state('fail1') == 'FAILED'
    assert node.get_task_state('after1') == 'SKIPPED'
    skipped_serialized = node.get_task_process('after1')
    assert skipped_serialized is None or deserialize_safe(skipped_serialized) is None

    serialized = node.get_task_process('fail1')
    assert serialized is not None
    excepted_process = deserialize_safe(serialized)
    assert isinstance(excepted_process, ProcessNode)
    assert excepted_process.process_state.value.upper() == 'EXCEPTED'
    # the recovered node is the same one CALL-linked into provenance
    assert [child.pk for child in node.called] == [excepted_process.pk]


def test_get_task_process_survives_repeated_call_link_labels() -> None:
    """A calcfunction inside a While zone reuses its call-link label every iteration.

    Regression test for aiidateam/aiida-workgraph#809: several CALL links can share
    a label when the same task runs more than once, so recovering the excepted node
    by exact-label lookup raises ``MultipleObjectsError``. The recovery must instead
    pick the newest matching node, keeping ``get_task_process`` pointed at the
    iteration that actually failed.
    """
    from aiida.orm import Int, ProcessNode
    from aiida_workgraph import While, WorkGraph, task
    from aiida_workgraph.orm.utils import deserialize_safe

    @task.calcfunction()
    def bump(x):
        if x.value >= 3:
            raise ValueError('deliberate failure on the third iteration')
        return x + 1

    @task()
    def smaller(x, y):
        return x < y

    with WorkGraph('test_while_repeated_label_process_reference') as wg:
        wg.ctx = {'n': Int(1)}
        cmp1 = wg.add_task(smaller, name='cmp1', x=wg.ctx.n, y=10)
        with While(cmp1.outputs.result, max_iterations=10) as zone:
            zone.add_task(bump, name='bump1', x=wg.ctx.n)
        wg.update_ctx({'n': wg.tasks.bump1.outputs.result})
        wg.run()

    node = wg.process
    serialized = node.get_task_process('bump1')
    assert serialized is not None
    excepted_process = deserialize_safe(serialized)
    assert isinstance(excepted_process, ProcessNode)
    assert excepted_process.process_state.value.upper() == 'EXCEPTED'
    # the recovered node is the newest 'bump1' CALL link, i.e. the failing
    # iteration, not an earlier successful iteration sharing the same label
    linked = node.base.links.get_outgoing(link_label_filter='bump1').all()
    assert excepted_process.pk == max(link.node.pk for link in linked)


def test_get_task_process_survives_missing_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    """A calcfunction run with ``metadata={'store_provenance': False}`` stores no node.

    Regression test for aiidateam/aiida-workgraph#809: with nothing CALL-linked into
    provenance, recovery has no node to find. It must degrade gracefully rather than
    raise its own lookup error, which would shadow the wrapped function's real
    exception in the engine's error log.
    """
    import logging

    from aiida.orm import Int
    from aiida_workgraph import WorkGraph, task
    from aiida_workgraph.orm.utils import deserialize_safe

    logged_errors = []
    original_error = logging.Logger.error

    def recording_error(self, msg, *args, **kwargs):
        logged_errors.append(str(msg))
        return original_error(self, msg, *args, **kwargs)

    monkeypatch.setattr(logging.Logger, 'error', recording_error)

    @task.calcfunction()
    def fails(x):
        raise ValueError('deliberate failure without provenance')

    wg = WorkGraph(name='test_noprovenance_process_reference')
    wg.add_task(fails, 'fail1', x=Int(1), metadata={'store_provenance': False})
    wg.run()

    node = wg.process
    assert node.exit_status == 302
    assert node.get_task_state('fail1') == 'FAILED'
    # nothing was stored or CALL-linked, so recovery has nothing to find
    assert node.called == []
    serialized = node.get_task_process('fail1')
    assert serialized is None or deserialize_safe(serialized) is None

    # the wrapped function's own exception reaches the engine's error log
    # unshadowed by a secondary error from the failed recovery attempt
    assert any('deliberate failure without provenance' in message for message in logged_errors)
