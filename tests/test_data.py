def test_PickledData():
    from aiida_workgraph.orm.pickled_data import PickledData

    class CustomData:
        def __init__(self, a):
            self.a = a

    data = CustomData(a=1)
    pickled_data = PickledData(data)
    pickled_data.store()
    assert pickled_data.value.a == 1


def test_PickledFunction():
    import numpy as np
    from aiida_workgraph.orm.pickled_function import PickledFunction

    data = PickledFunction.inspect_function(np.add)
    assert data == {
        "name": "add",
        "source_code": "",
        "source_code_without_decorator": "",
        "import_statements": "",
    }
