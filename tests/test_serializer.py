import aiida

aiida.load_profile()


def test_python_job():
    """Test a simple python node."""
    from aiida_workgraph.orm import GeneralData, general_serializer

    inputs = {"a": 1, "b": 2.0, "c": set()}
    new_inputs = general_serializer(inputs)
    assert isinstance(new_inputs["a"], aiida.orm.Int)
    assert isinstance(new_inputs["b"], aiida.orm.Float)
    assert isinstance(new_inputs["c"], GeneralData)
