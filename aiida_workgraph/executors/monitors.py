async def monitor(function, interval, *args, **kwargs):
    import asyncio

    while True:
        result = function(*args, **kwargs)
        if result:
            break
        await asyncio.sleep(interval)


def file_monitor(filename):
    """Check if the file exists."""
    import os

    return os.path.exists(filename)


def time_monitor(time):
    """Return True if the current time is greater than the given time."""
    import datetime

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
    state = node.base.extras.get(f"_task_state_{task_name}", None)
    print(f"Task state: {state}")
    if state in ["FINISHED", "FAILED", "SKIPPED"]:
        return True
    else:
        return False
