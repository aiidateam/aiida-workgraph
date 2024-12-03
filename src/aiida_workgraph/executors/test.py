from typing import Union, Dict
import time
from aiida.engine import calcfunction
from aiida.orm import Int, Float
from aiida_workgraph import task


@calcfunction
def add(
    x: Union[Int, Float], y: Union[Int, Float], t: Union[Int, Float] = 1.0
) -> Dict[str, Union[Int, Float]]:
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y}


@calcfunction
def sum_diff(
    x: Union[Int, Float], y: Union[Int, Float], t: Union[Int, Float] = 1.0
) -> Dict[str, Union[Int, Float]]:
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y, "diff": x - y}


@task.pythonjob()
def add_pythonjob(x: int, y: int) -> int:
    return x + y
