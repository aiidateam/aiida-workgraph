import aiida


def test_python_job():
    """Test a simple python node."""
    from aiida_workgraph.orm import GeneralData, serialize_to_aiida_nodes

    inputs = {"a": 1, "b": 2.0, "c": set()}
    new_inputs = serialize_to_aiida_nodes(inputs)
    assert isinstance(new_inputs["a"], aiida.orm.Int)
    assert isinstance(new_inputs["b"], aiida.orm.Float)
    assert isinstance(new_inputs["c"], GeneralData)
