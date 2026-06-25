import logging
import pytest
import time

from aiida_workgraph.engine.task_actions import TaskActionManager
from aiida_workgraph.enums import TaskAction, TaskState


class _RecordingStateManager:
    """Minimal stand-in for ``TaskStateManager`` that records dispatched calls.

    A real one needs a live engine/process node, so the collaborator is faked to
    keep the dispatch test fast and focused on ``apply_task_actions``.
    """

    def __init__(self):
        self.calls = []

    def reset_task(self, name):
        self.calls.append(('reset', name))

    def set_task_runtime_info(self, name, key, value):
        self.calls.append(('set', name, key, value))


class _RecordingProcess:
    def __init__(self):
        self.reports = []

    def report(self, message):
        self.reports.append(message)


def _make_action_manager():
    state_manager = _RecordingStateManager()
    process = _RecordingProcess()
    manager = TaskActionManager(state_manager, logging.getLogger(__name__), process)
    return manager, state_manager, process


@pytest.mark.parametrize(
    'raw_action, expected_call',
    [
        pytest.param('pause', ('set', 'add1', 'action', TaskAction.PAUSE), id='pause-lowercase'),
        pytest.param('PAUSE', ('set', 'add1', 'action', TaskAction.PAUSE), id='pause-uppercase'),
        pytest.param('PaUsE', ('set', 'add1', 'action', TaskAction.PAUSE), id='pause-mixedcase'),
        pytest.param('reset', ('reset', 'add1'), id='reset-lowercase'),
        pytest.param('ReSeT', ('reset', 'add1'), id='reset-mixedcase'),
        pytest.param('skip', ('set', 'add1', 'state', TaskState.SKIPPED), id='skip-lowercase'),
    ],
)
def test_apply_task_actions_is_case_insensitive(raw_action, expected_call):
    """Whatever the case of the incoming action, it dispatches the same way."""
    manager, state_manager, _ = _make_action_manager()
    manager.apply_task_actions({'action': raw_action, 'tasks': ['add1']})
    assert expected_call in state_manager.calls


def test_apply_task_actions_raises_on_unknown_action():
    """A typo'd action fails loudly rather than being silently swallowed."""
    manager, state_manager, _ = _make_action_manager()
    with pytest.raises(ValueError):
        manager.apply_task_actions({'action': 'frobnicate', 'tasks': ['add1']})
    assert state_manager.calls == []


@pytest.mark.skip(reason='PAUSED state is wrong for the moment.')
def test_pause_play_workgraph(wg_engine):
    wg = wg_engine
    wg.name = 'test_pause_play_workgraph'
    wg.submit()
    time.sleep(5)
    wg.pause()
    wg.update()
    assert wg.process.process_state.value.upper() == 'PAUSED'


# @pytest.mark.skip(reason='pause task is not stable for the moment.')
@pytest.mark.usefixtures('started_daemon_client')
def test_pause_play_task(wg_calcjob):
    wg = wg_calcjob
    wg.name = 'test_pause_play_task'
    # pause add1 before submit
    wg.pause_tasks(['add1'])
    wg.submit()
    # wait for the workgraph to launch add1
    wg.wait(tasks={'add1': ['CREATED']}, timeout=60, interval=5)
    assert wg.tasks.add1.node.process_state.value.upper() == 'CREATED'
    assert wg.tasks.add1.node.process_status == 'Paused through WorkGraph'
    # pause add2 after submit
    wg.pause_tasks(['add2'])
    time.sleep(5)
    wg.play_tasks(['add1'])
    # wait for the workgraph to launch add2
    wg.wait(tasks={'add2': ['CREATED']}, timeout=60, interval=5)
    assert wg.tasks.add2.node.process_state.value.upper() == 'CREATED'
    assert wg.tasks.add2.node.process_status == 'Paused through WorkGraph'
    # I disabled the following lines because the test is not stable
    # Seems the daemon is not responding to the play signal
    wg.play_tasks(['add2'])
    wg.wait(interval=5)
    assert wg.tasks.add2.outputs.sum.value == 9


def test_pause_play_error_handler(wg_calcjob, create_process_node):
    wg = wg_calcjob
    wg.name = 'test_pause_play_error_handler'
    wg.process = create_process_node(state='finished', exit_status=0)
    try:
        wg.pause_tasks(['add1'])
    except Exception as e:
        assert 'WorkGraph is finished. Cannot pause tasks.' in str(e)

    try:
        wg.play_tasks(['add1'])
    except Exception as e:
        assert 'WorkGraph is finished. Cannot play tasks.' in str(e)

    try:
        wg.kill_tasks(['add2'])
    except Exception as e:
        assert 'WorkGraph is finished. Cannot kill tasks.' in str(e)
