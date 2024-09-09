import pytest
from aiida_workgraph import WorkGraph, task
from aiida import orm


@pytest.mark.parametrize(
    "data_type, socket_type",
    (
        (int, "workgraph.int"),
        (float, "workgraph.float"),
        (bool, "workgraph.bool"),
        (str, "workgraph.string"),
        (orm.Int, "workgraph.aiida_int"),
        (orm.Float, "workgraph.aiida_float"),
        (orm.Str, "workgraph.aiida_string"),
        (orm.Bool, "workgraph.aiida_bool"),
    ),
)
def test_type_mapping(data_type, socket_type) -> None:
    """Test the mapping of data types to socket types."""

    @task()
    def add(x: data_type):
        pass

    assert add.task().inputs["x"].identifier == socket_type


def test_socket(decorated_multiply) -> None:
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
    # wg.submit(wait=True)
    wg.run()
    assert wg.state.upper() == "FINISHED"
    assert multiply1.outputs["result"].value == 6.0


def test_AiiDA_socket():
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


@pytest.mark.parametrize(
    "data_type, socket_value, node_value",
(
    (None, None, None),
    # Check that SocketAny works for int node, without providing type hint
    (None, 1, 1),
    (int, 1, 1),
    (float, 1., 1.),
    (bool, True, True),
    (str, "abc", "abc"),
    # TODO: Wanted to instantiate the AiiDA ORM classes here, but that raises an
    # TODO: aiida.common.exceptions.ConfigurationError, due to profile not being loaded
    # TODO: which also isn't resolved by: `@pytest.mark.usefixtures("aiida_profile")`
    (orm.Int, 1, 1),
    (orm.Float, 1., 1.),
    (orm.Str, 'abc', 'abc'),
    (orm.Bool, True, True),
))
def test_node_value(data_type, socket_value, node_value):

    wg = WorkGraph()

    def my_task(x: data_type):
        pass

    my_task1 = wg.add_task(my_task, name="my_task", x=socket_value)
    socket = my_task1.inputs['x']

    # Private attribute is undefined (and shouldn't be called anyway)
    with pytest.raises(AttributeError):
        socket._node_value

    # This should call the `get_node_value` method
    socket_node_value = socket.node_value

    assert isinstance(socket_node_value, type(node_value))
    assert socket_node_value == node_value

    # Now the private attribute should be set, such that the `get_node_value` method doesn't have to be called again
    assert isinstance(socket_node_value, type(socket._node_value))
    assert socket_node_value == socket._node_value
