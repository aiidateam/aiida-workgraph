import pytest
from aiida_workgraph import WorkGraph, task
from aiida import orm
from typing import Any


# @pytest.mark.usefixtures("aiida_profile")
@pytest.mark.parametrize(
    "data_type, data, identifier",
    (
        (Any, 1, "workgraph.any"),
        (int, 1, "workgraph.int"),
        (float, 2.0, "workgraph.float"),
        (bool, True, "workgraph.bool"),
        (str, "abc", "workgraph.string"),
        (str, "{{variable}}", "workgraph.string"),
        (orm.Int, 1, "workgraph.aiida_int"),
        (orm.Float, 2.0, "workgraph.aiida_float"),
        (orm.Str, "abc", "workgraph.aiida_string"),
        (orm.Bool, True, "workgraph.aiida_bool"),
        (orm.Bool, "{{variable}}", "workgraph.aiida_bool"),
    ),
)
def test_type_mapping(data_type, data, identifier) -> None:
    """Test the mapping of data types to socket types."""

    @task()
    def add(x: data_type):
        pass

    assert add.task().inputs["x"].identifier == identifier
    assert add.task().inputs["x"].property.identifier == identifier
    add_task = add.task()
    add_task.set({"x": data})
    # test set data from context
    add_task.set({"x": "{{variable}}"})


def test_aiida_data_socket() -> None:
    """Test the mapping of data types to socket types."""
    from aiida import orm, load_profile

    load_profile()

    datas = [(orm.StructureData, orm.StructureData(), "workgraph.aiida_structuredata")]
    for data_type, data, identifier in datas:

        @task()
        def add(x: data_type):
            pass

        assert add.task().inputs["x"].identifier == identifier
        assert add.task().inputs["x"].property.identifier == identifier
        add_task = add.task()
        add_task.set({"x": data})
        # test set data from context
        add_task.set({"x": "{{variable}}"})


@pytest.mark.parametrize(
    "data_type, data",
    (
        (int, 1.0),
        (float, "a"),
        (bool, "a"),
        (str, [1, 2, 3]),
        (orm.Int, 1.0),
        (orm.Float, "a"),
        (orm.Str, [1, 2, 3]),
        (orm.Bool, "a"),
        (orm.StructureData, 1),
    ),
)
def test_socket_validate(data_type, data) -> None:
    @task()
    def add(x: data_type):
        """"""

    add_task = add.task()
    # Test setting a value that should raise an exception
    with pytest.raises(Exception) as excinfo:
        add_task.set({"x": data})

    assert "Expected value of type" in str(excinfo.value)


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


def test_kwargs() -> None:
    """Test the kwargs of a task."""

    @task()
    def test(a, b=1, **kwargs):
        return {"sum": a + b, "product": a * b}

    test1 = test.node()
    assert test1.inputs["kwargs"].link_limit == 1e6
    assert test1.inputs["kwargs"].identifier == "workgraph.namespace"
    assert test1.inputs["kwargs"].property.value is None


@pytest.mark.parametrize(
    "data_type, socket_value, node_value",
    (
        (None, None, None),
        # Check that SocketAny works for int node, without providing type hint
        (None, 1, 1),
        (int, 1, 1),
        (float, 1.0, 1.0),
        (bool, True, True),
        (str, "abc", "abc"),
        (orm.Int, 1, 1),
        (orm.Float, 1.0, 1.0),
        (orm.Str, "abc", "abc"),
        (orm.Bool, True, True),
    ),
)
def test_node_value(data_type, socket_value, node_value):

    wg = WorkGraph()

    def my_task(x: data_type):
        pass

    my_task1 = wg.add_task(my_task, name="my_task", x=socket_value)
    socket = my_task1.inputs["x"]

    socket_node_value = socket.get_node_value()
    assert isinstance(socket_node_value, type(node_value))
    assert socket_node_value == node_value

    # Check that property also returns the correct results
    socket_node_value = socket.node_value
    assert isinstance(socket_node_value, type(node_value))
    assert socket_node_value == node_value
