from typing import Union, Dict
import time
from aiida.engine import calcfunction
from aiida.orm import Int, Float


@calcfunction
def add(
    x: Union[Int, Float], y: Union[Int, Float], t: Union[Int, Float] = Float(1.0)
) -> Dict[str, Union[Int, Float]]:
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y}


@calcfunction
def greater(
    x: Union[Int, Float], y: Union[Int, Float], t: Union[Int, Float] = Float(1.0)
) -> Dict[str, bool]:
    """Compare node."""
    time.sleep(t.value)
    return {"result": x > y}


@calcfunction
def sum_diff(
    x: Union[Int, Float], y: Union[Int, Float], t: Union[Int, Float] = Float(1.0)
) -> Dict[str, Union[Int, Float]]:
    """Add node."""
    time.sleep(t.value)
    return {"sum": x + y, "diff": x - y}
