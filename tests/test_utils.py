import pytest
from aiida import orm
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


def test_get_parent_workgraphs():
    from aiida.common.links import LinkType
    from aiida_workgraph.utils import get_parent_workgraphs

    wn1 = orm.WorkflowNode()
    wn2 = orm.WorkflowNode()
    wn3 = orm.WorkflowNode()
    wn3.base.links.add_incoming(wn2, link_type=LinkType.CALL_WORK, link_label="link")
    wn2.base.links.add_incoming(wn1, link_type=LinkType.CALL_WORK, link_label="link")
    wn1.store()
    wn2.store()
    wn3.store()

    parent_workgraphs = get_parent_workgraphs(wn3.pk)
    assert len(parent_workgraphs) == 3


def test_generate_node_graph():
    from aiida_workgraph.utils import generate_node_graph
    from IPython.display import IFrame
    import os

    wn1 = orm.WorkflowNode()
    wn1.store()

    graph = generate_node_graph(wn1.pk)
    assert isinstance(graph, IFrame)
    # check file html/node_graph_{pk}.html is created
    assert os.path.isfile(f"html/node_graph_{wn1.pk}.html")
