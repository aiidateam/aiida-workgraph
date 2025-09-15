from aiida_workgraph import task, namespace
from typing import Annotated
import pytest
import re


@task
def add(x, y):
    return x + y


def test_validate_required_inputs():
    @task.graph()
    def my_graph(a, b: Annotated[dict, namespace(x=int, y=int)]):
        add(a, b["x"])
        add(a)

    with pytest.raises(
        ValueError,
        match=re.escape("Missing required inputs: ['graph_inputs.b.y', 'add1.y']"),
    ):
        my_graph.run(a=1, b={"x": 1})
