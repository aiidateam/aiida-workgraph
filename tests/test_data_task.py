import pytest
from aiida_workgraph import WorkGraph


@pytest.mark.parametrize(
    "identifier, data",
    (
        ("workgraph.aiida_int", 1),
        ("workgraph.aiida_float", 2.0),
        ("workgraph.aiida_string", "abc"),
    ),
)
def test_data_task(identifier, data) -> None:
    """Test a normal task."""

    wg = WorkGraph("test_normal_task")
    task1 = wg.add_task(identifier, name="task1", value=data)
    wg.run()
    assert task1.outputs["result"].value == data


def test_data_dict_task():
    """Test a normal task."""

    wg = WorkGraph("test_data_dict_task")
    task1 = wg.add_task("workgraph.aiida_dict", name="task1", value={"a": 1})
    wg.run()
    assert task1.outputs["result"].value == {"a": 1}


def test_data_list_task():
    """Test a normal task."""

    wg = WorkGraph("test_data_list_task")
    task1 = wg.add_task("workgraph.aiida_list", name="task1", value=[1, 2, 3])
    wg.run()
    assert task1.outputs["result"].value == [1, 2, 3]
