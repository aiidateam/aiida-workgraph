import time
from aiida.engine import calcfunction


@calcfunction
def add(x, y, t=1):
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y}


@calcfunction
def greater(x, y, t=1):
    """Compare node."""
    time.sleep(t.value)
    return {"result": x > y}


@calcfunction
def sum_diff(x, y, t=1.0):
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y, "diff": x - y}


if __name__ == "__main__":
    pass
