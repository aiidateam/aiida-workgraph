import pytest
from aiida_workgraph import WorkGraph


def test_socket(decorated_multiply) -> None:
    from aiida_workgraph import task

    @task(
        outputs=[{"identifier": float, "name": "result"}],
    )
    def add(x: int, y: float):
        result = x + y
        return result

    #
    wg = WorkGraph()
    add1 = wg.add_task(add, name="add1")
    multiply1 = wg.add_task(
        decorated_multiply, name="multiply1", x=add1.outputs["result"], y=2
    )
    # Test setting a value that should raise an exception
    with pytest.raises(Exception) as excinfo:
        add1.set({"x": 1.0, "y": 2})  # Direct integers instead of orm.Int or orm.Float

    assert "is not" in str(excinfo.value) and int.__name__ in str(excinfo.value)
    # This should be successful
    add1.set({"x": 1, "y": 2.0})
    wg.submit(wait=True)
    assert wg.state.upper() == "FINISHED"
    assert multiply1.outputs["result"].value == 6.0


def test_AiiDA_socket():
    from aiida_workgraph import task
    from aiida import orm

    @task.calcfunction(
        outputs=[{"identifier": orm.Float, "name": "result"}],
    )
    def add(x: int, y: float):
        result = x + y
        return result

    wg = WorkGraph()
    add1 = wg.add_task(add, name="add1")
    # Test setting a value that should raise an exception
    with pytest.raises(Exception) as excinfo:
        add1.set(
            {"x": orm.Float(1.0), "y": 2}
        )  # Direct integers instead of orm.Int or orm.Float

    assert "is not" in str(excinfo.value)
    # This should be successful
    add1.set({"x": orm.Int(1), "y": orm.Float(2.0)})
    wg.run()
    assert wg.state.upper() == "FINISHED"
    assert wg.tasks["add1"].outputs["result"].value == 3.0


@pytest.mark.usefixtures("started_daemon_client")
def test_numpy_array(decorated_normal_add):
    """Test data type with numpy array."""
    import numpy as np

    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    wg = WorkGraph()
    wg.add_task(decorated_normal_add, name="add1", x=x, y=y)
    wg.submit(wait=True)
    # wg.run()
    assert wg.state.upper() == "FINISHED"
    # assert (wg.tasks["add1"].outputs["result"].value == np.array([5, 7, 9])).all()
