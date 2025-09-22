import datetime
import typing as t

from aiida_workgraph import task


@task.monitor
def monitor_file(filepath: str):
    """Return `True` when the file is detected."""
    import os

    return os.path.exists(filepath)


@task.monitor
def monitor_time(time: t.Union[str, datetime.datetime]):
    """Return `True` when the given moment in time has passed.

    If given as a string, `time` should be in ISO format (e.g., '2025-07-30T12:00:00').
    """

    if isinstance(time, str):
        try:
            time = datetime.datetime.fromisoformat(time)
        except ValueError as err:
            raise ValueError(f'Invalid time format: {time}. Expected ISO format.') from err

    return datetime.datetime.now() > time


@task.monitor
def monitor_task(task_name: str, workgraph_pk: int = None, workgraph_name: str = None):
    """Return `True` if the task in the WorkGraph is completed."""
    from aiida import orm

    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    if workgraph_pk:
        try:
            node = orm.load_node(workgraph_pk)
        except Exception:
            return False
    else:
        builder = orm.QueryBuilder()
        builder.append(
            WorkGraphEngine,
            filters={'attributes.process_label': {'==': f'WorkGraph<{workgraph_name}>'}},
            tag='process',
        )
        if builder.count() == 0:
            return False
    print('Found workgraph')
    node = builder.first()[0]
    state = node.task_states.get(task_name, '')
    print(f'Task state: {state}')
    return state in ['FINISHED', 'FAILED', 'SKIPPED']
