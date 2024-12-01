import pytest

from aiida_workgraph.utils import validate_task_inout


def test_validate_task_inout_str_list():
    """Test validation with a list of strings."""
    input_list = ["task1", "task2"]
    result = validate_task_inout(input_list, "input")
    assert result == [{"name": "task1"}, {"name": "task2"}]


def test_validate_task_inout_dict_list():
    """Test validation with a list of dictionaries."""
    input_list = [{"name": "task1"}, {"name": "task2"}]
    result = validate_task_inout(input_list, "input")
    assert result == input_list


@pytest.mark.parametrize(
    "input_list, list_type, expected_error",
    [
        # Mixed types error cases
        (
            ["task1", {"name": "task2"}],
            "input",
            "Provide either a list of `str` or `dict` as `input`, not mixed types.",
        ),
        (
            [{"name": "task1"}, "task2"],
            "output",
            "Provide either a list of `str` or `dict` as `output`, not mixed types.",
        ),
        # Empty list cases
        ([], "input", None),
        ([], "output", None),
    ],
)
def test_validate_task_inout_mixed_types(input_list, list_type, expected_error):
    """Test error handling for mixed type lists."""
    if expected_error:
        with pytest.raises(TypeError) as excinfo:
            validate_task_inout(input_list, list_type)
        assert str(excinfo.value) == expected_error
    else:
        # For empty lists, no error should be raised
        result = validate_task_inout(input_list, list_type)
        assert result == []


@pytest.mark.parametrize(
    "input_list, list_type",
    [
        # Invalid type cases
        ([1, 2, 3], "input"),
        ([None, None], "output"),
        ([True, False], "input"),
        (["task", 123], "output"),
    ],
)
def test_validate_task_inout_invalid_types(input_list, list_type):
    """Test error handling for completely invalid type lists."""
    with pytest.raises(TypeError) as excinfo:
        validate_task_inout(input_list, list_type)
    assert "Provide either a list of" in str(excinfo.value)


def test_validate_task_inout_dict_with_extra_keys():
    """Test validation with dictionaries having extra keys."""
    input_list = [
        {"name": "task1", "description": "first task"},
        {"name": "task2", "priority": "high"},
    ]
    result = validate_task_inout(input_list, "input")
    assert result == input_list


def test_get_or_create_code(fixture_localhost):
    from aiida_workgraph.utils import get_or_create_code
    from aiida.orm import Code

    # create a new code
    code1 = get_or_create_code(
        computer="localhost",
        code_label="test_code",
        code_path="/bin/bash",
        prepend_text='echo "Hello, World!"',
    )
    assert isinstance(code1, Code)
    # use already created code
    code2 = get_or_create_code(
        computer="localhost",
        code_label="test_code",
        code_path="/bin/bash",
        prepend_text='echo "Hello, World!"',
    )
    assert code1.uuid == code2.uuid
