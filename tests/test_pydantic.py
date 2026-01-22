import pytest

from aiida_workgraph import task

try:
    from pydantic import BaseModel as _BaseModel
except Exception:
    _BaseModel = None

if _BaseModel is not None:

    class PydanticInputs(_BaseModel):
        x: int
        y: int

    class PydanticOutputs(_BaseModel):
        sum: int
        product: int

else:
    PydanticInputs = None
    PydanticOutputs = None


if _BaseModel is not None:

    @task
    def add_multiply_pydantic(data: 'PydanticInputs') -> 'PydanticOutputs':
        return PydanticOutputs(sum=data.x + data.y, product=data.x * data.y)
else:
    add_multiply_pydantic = None


def test_workgraph_pydantic_inputs_outputs():
    pytest.importorskip('pydantic')

    @task.graph
    def add_graph(data: 'PydanticInputs') -> 'PydanticOutputs':
        return add_multiply_pydantic(data=data)

    result, wg = add_graph.run_get_graph(data=PydanticInputs(x=2, y=3))

    assert result['sum'] == 5
    assert result['product'] == 6
    assert wg.outputs.sum.value == 5
    assert wg.outputs.product.value == 6
