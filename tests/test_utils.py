import pytest

from aiida_workgraph.utils import validate_task_inout


def test_validate_task_inout_empty_list():
    """Test validation with a list of strings."""
    input_list = []
    result = validate_task_inout(input_list, "inputs")
    assert result == []


def test_validate_task_inout_str_list():
    """Test validation with a list of strings."""
    input_list = ["task1", "task2"]
    result = validate_task_inout(input_list, "inputs")
    assert result == [{"name": "task1"}, {"name": "task2"}]


def test_validate_task_inout_dict_list():
    """Test validation with a list of dictionaries."""
    input_list = [{"name": "task1"}, {"name": "task2"}]
    result = validate_task_inout(input_list, "inputs")
    assert result == input_list


def test_validate_task_inout_mixed_list():
    """Test validation with a list of dictionaries."""
    input_list = ["task1", {"name": "task2"}]
    result = validate_task_inout(input_list, "inputs")
    assert result == [{"name": "task1"}, {"name": "task2"}]


@pytest.mark.parametrize(
    "input_list, list_type",
    [
        # Invalid type cases
        ([1, 2, 3], "inputs"),
        ([None, None], "outputs"),
        ([True, False], "inputs"),
        (["task", 123], "outputs"),
    ],
)
def test_validate_task_inout_invalid_types(input_list, list_type):
    """Test error handling for completely invalid type lists."""
    with pytest.raises(TypeError) as excinfo:
        validate_task_inout(input_list, list_type)
    assert "Wrong type provided" in str(excinfo.value)


def test_validate_task_inout_dict_with_extra_keys():
    """Test validation with dictionaries having extra keys."""
    input_list = [
        {"name": "task1", "description": "first task"},
        {"name": "task2", "priority": "high"},
    ]
    result = validate_task_inout(input_list, "inputs")
    assert result == input_list
