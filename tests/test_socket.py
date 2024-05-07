import pytest
from aiida_workgraph import WorkGraph
import aiida

aiida.load_profile()


def test_socket(decorated_multiply) -> None:
    from aiida_workgraph import node

    @node(
        inputs=[[int, "x"], [float, "y"]],
        outputs=[[float, "result"]],
    )
    def add(x, y):
        result = x + y
        return result

    #
    wg = WorkGraph()
    add1 = wg.nodes.new(add, name="add1")
    multiply1 = wg.nodes.new(
        decorated_multiply, name="multiply1", x=add1.outputs["result"], y=2
    )
    # Test setting a value that should raise an exception
    with pytest.raises(Exception) as excinfo:
        add1.set({"x": 1.0, "y": 2})  # Direct integers instead of orm.Int or orm.Float

    assert "is not an {}".format(int.__name__) in str(excinfo.value)
    # This should be successful
    add1.set({"x": 1, "y": 2.0})
    wg.submit(wait=True)
    assert wg.state.upper() == "FINISHED"
    assert multiply1.outputs["result"].value == 6.0


def test_AiiDA_socket():
    from aiida_workgraph import node
    from aiida import orm

    @node.calcfunction(
        inputs=[[orm.Int, "x"], [orm.Float, "y"]],
        outputs=[[orm.Float, "result"]],
    )
    def add(x, y):
        result = x + y
        return result

    wg = WorkGraph()
    add1 = wg.nodes.new(add, name="add1")
    # Test setting a value that should raise an exception
    with pytest.raises(Exception) as excinfo:
        add1.set(
            {"x": orm.Float(1.0), "y": 2}
        )  # Direct integers instead of orm.Int or orm.Float

    assert "is not an {}".format(orm.Int.__name__) in str(excinfo.value)
    # This should be successful
    add1.set({"x": orm.Int(1), "y": orm.Float(2.0)})
    wg.run()
    assert wg.state.upper() == "FINISHED"
    assert wg.nodes["add1"].outputs["result"].value == 3.0


def test_numpy_array(decorated_normal_add):
    """Test data type with numpy array."""
    import numpy as np

    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    wg = WorkGraph()
    wg.nodes.new(decorated_normal_add, name="add1", x=x, y=y)
    wg.submit(wait=True)
    # wg.run()
    assert wg.state.upper() == "FINISHED"
    # assert (wg.nodes["add1"].outputs["result"].value == np.array([5, 7, 9])).all()
