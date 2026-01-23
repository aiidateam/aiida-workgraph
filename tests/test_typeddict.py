import pytest
from typing import TypedDict

from aiida_workgraph import task

try:
    from typing import Unpack
except ImportError:
    try:
        from typing_extensions import Unpack
    except Exception:
        Unpack = None


class XYIn(TypedDict):
    x: int
    y: int


class AddMultiplyOut(TypedDict):
    sum: int
    product: int


@task(inputs=XYIn)
def add_multiply_typed_dict(**kwargs) -> AddMultiplyOut:
    return {'sum': kwargs['x'] + kwargs['y'], 'product': kwargs['x'] * kwargs['y']}


if Unpack is not None:

    @task
    def add_multiply_typed_dict_unpack(**kwargs: 'Unpack[XYIn]') -> AddMultiplyOut:
        return {'sum': kwargs['x'] + kwargs['y'], 'product': kwargs['x'] * kwargs['y']}
else:
    add_multiply_typed_dict_unpack = None


def test_workgraph_typeddict_inputs_outputs():
    @task.graph
    def add_graph(x: int, y: int) -> AddMultiplyOut:
        return add_multiply_typed_dict(x=x, y=y)

    result, wg = add_graph.run_get_graph(x=2, y=3)

    assert result['sum'] == 5
    assert result['product'] == 6
    assert wg.outputs.sum.value == 5
    assert wg.outputs.product.value == 6


def test_workgraph_typeddict_unpack_kwargs():
    if Unpack is None or add_multiply_typed_dict_unpack is None:
        pytest.skip('Unpack not available')

    @task.graph
    def add_graph(x: int, y: int) -> AddMultiplyOut:
        return add_multiply_typed_dict_unpack(x=x, y=y)

    result, wg = add_graph.run_get_graph(x=2, y=3)

    assert result['sum'] == 5
    assert result['product'] == 6
    assert wg.outputs.sum.value == 5
    assert wg.outputs.product.value == 6
