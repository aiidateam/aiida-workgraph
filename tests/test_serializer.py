"""
Test serialize_input_values_recursively and deserialize_input_values_recursively
"""
from aiida_workgraph import WorkGraph, task


@task.graph_builder()
def sub_workflow(func):

    wg = WorkGraph("sub_workflow")
    wg.add_task(func)
    return wg


def test_func_as_input(capsys):
    from aiida_workgraph.executors.test import add

    wg = WorkGraph("test_func_as_input")
    wg.add_task(sub_workflow, func=add, name="sub_workflow")
    wg.save()

    # load and capture stdout
    loaded_wg = WorkGraph.load(wg.pk)
    captured = capsys.readouterr()

    assert "Info: could not deserialize input" in captured.out
    assert "sub_workflow" in loaded_wg.tasks
