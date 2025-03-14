def test_PickledData():
    from aiida_workgraph.orm.pickled_data import PickledData

    class CustomData:
        def __init__(self, a):
            self.a = a

    data = CustomData(a=1)
    pickled_data = PickledData(data)
    pickled_data.store()
    assert pickled_data.value.a == 1
