import time
from aiida.orm import Int
from aiida_workgraph import task
from aiida_workgraph.socket_spec import namespace


@task.calcfunction
def add(x: Int = 0, y: Int = 0, t: Int = 1) -> Int:
    """Add node."""
    time.sleep(t.value)
    return x + y


@task.calcfunction
def sum_diff(x: Int = 0, y: Int = 0, t: Int = 1) -> namespace(sum=Int, diff=Int):
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y, "diff": x - y}


@task.pythonjob()
def add_pythonjob(x: int, y: int) -> int:
    return x + y
