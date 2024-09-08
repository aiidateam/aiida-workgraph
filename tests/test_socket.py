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
        # (orm.StructureData, orm.StructureData(), "workgraph.aiida_structuredata"),
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
