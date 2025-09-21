from __future__ import annotations

from aiida.manage import get_manager
from aiida import orm
from aiida.engine.processes import control


def create_task_action(
    pk: int,
    tasks: list,
    action: str = 'pause',
):
    """Send task action to Process."""

    controller = get_manager().get_process_controller()
    message = {'intent': 'custom', 'catalog': 'task', 'action': action, 'tasks': tasks}
    controller._communicator.rpc_send(pk, message)


def get_task_runtime_info(node, name: str, key: str) -> str:
    """Get task state info from attributes."""
    from aiida_workgraph.orm.utils import deserialize_safe

    if key == 'process':
        value = deserialize_safe(node.task_processes.get(name, ''))
    elif key == 'state':
        value = node.task_states.get(name, '')
    elif key == 'action':
        value = node.task_actions.get(name, '')
    return value


def pause_tasks(pk: int, tasks: list[str], timeout: int = 5):
    """Pause task."""
    node = orm.load_node(pk)
    if node.is_finished:
        message = 'WorkGraph is finished. Cannot pause tasks.'
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        'CREATED',
        'RUNNING',
        'WAITING',
        'PAUSED',
    ]:
        for name in tasks:
            if get_task_runtime_info(node, name, 'state') == 'PLANNED':
                create_task_action(pk, tasks, action='pause')
            elif get_task_runtime_info(node, name, 'state') == 'RUNNING':
                try:
                    control.pause_processes(
                        [get_task_runtime_info(node, name, 'process')],
                        all_entries=None,
                        timeout=timeout,
                    )
                except Exception as e:
                    print(f'Pause task {name} failed: {e}')
    return True, ''


def play_tasks(pk: int, tasks: list, timeout: int = 5):
    node = orm.load_node(pk)
    if node.is_finished:
        message = 'WorkGraph is finished. Cannot kill tasks.'
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        'CREATED',
        'RUNNING',
        'WAITING',
        'PAUSED',
    ]:
        for name in tasks:
            state = get_task_runtime_info(node, name, 'state')
            if state == 'PLANNED':
                create_task_action(pk, tasks, action='play')
                break
            process = get_task_runtime_info(node, name, 'process')
            if process.is_finished:
                raise ValueError(f'Task {name} is already finished.')
            elif process.process_state.value.upper() in ['CREATED', 'WAITING']:
                try:
                    control.play_processes(
                        [process],
                        all_entries=None,
                        timeout=timeout,
                    )
                except Exception as e:
                    print(f'Play task {name} failed: {e}')
    return True, ''


def kill_tasks(pk: int, tasks: list, timeout: int = 5):
    node = orm.load_node(pk)
    if node.is_finished:
        message = 'WorkGraph is finished. Cannot kill tasks.'
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        'CREATED',
        'RUNNING',
        'WAITING',
        'PAUSED',
    ]:
        for name in tasks:
            state = get_task_runtime_info(node, name, 'state')
            process = get_task_runtime_info(node, name, 'process')
            if state == 'PLANNED':
                create_task_action(pk, tasks, action='skip')
            elif state in [
                'CREATED',
                'RUNNING',
                'WAITING',
                'PAUSED',
            ]:
                if process is None:
                    print(f'Task {name} is not a AiiDA process.')
                    create_task_action(pk, tasks, action='kill')
                else:
                    try:
                        control.kill_processes(
                            [process],
                            all_entries=None,
                            timeout=timeout,
                        )
                    except Exception as e:
                        print(f'Kill task {name} failed: {e}')
    return True, ''


def reset_tasks(pk: int, tasks: list) -> None:
    """Reset tasks
    Args:
        tasks (list): a list of task names.
    """
    node = orm.load_node(pk)
    if node.is_finished:
        message = 'WorkGraph is finished. Cannot kill tasks.'
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        'CREATED',
        'RUNNING',
        'WAITING',
        'PAUSED',
    ]:
        for name in tasks:
            create_task_action(pk, tasks, action='reset')

    return True, ''
