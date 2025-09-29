from aiida_workgraph import WorkGraph, task


@task.graph()
def sub_workflow(func):
    func()


def test_func_as_input(capsys):
    from aiida_workgraph.executors.test import add

    wg = WorkGraph('test_func_as_input')
    wg.add_task(sub_workflow, func=add, name='sub_workflow')
    wg.save()

    # load and capture stdout
    loaded_wg = WorkGraph.load(wg.pk)
    captured = capsys.readouterr()

    assert 'Info: could not deserialize input' in captured.out
    assert 'sub_workflow' in loaded_wg.tasks
