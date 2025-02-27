async def monitor(function, interval, timeout, *args, **kwargs):
    """Monitor the function until it returns True or the timeout is reached."""
    import asyncio
    import time

    start_time = time.time()
    while True:
        result = function(*args, **kwargs)
        if result:
            break
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Timeout reached for monitor function {function}")
        await asyncio.sleep(interval)


def file_monitor(filepath: str):
    """Check if the file exists."""
    import os

    return os.path.exists(filepath)


def time_monitor(time: str):
    """Return True if the current time is greater than the given time."""
    import datetime

    # load the time string
    time = datetime.datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f")

    return datetime.datetime.now() > time


def task_monitor(task_name: str, workgraph_pk: int = None, workgraph_name: str = None):
    """Return True if the task in the WorkGraph is completed."""
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
            filters={
                "attributes.process_label": {"==": f"WorkGraph<{workgraph_name}>"}
            },
            tag="process",
        )
        if builder.count() == 0:
            return False
    print("Found workgraph")
    node = builder.first()[0]
    state = node.task_states.get(task_name, "")
    print(f"Task state: {state}")
    if state in ["FINISHED", "FAILED", "SKIPPED"]:
        return True
    else:
        return False
