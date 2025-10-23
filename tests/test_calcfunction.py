from aiida_workgraph import WorkGraph, task
from aiida import orm


def test_run(wg_task: WorkGraph) -> None:
    """Run simple calcfunction."""
    wg = wg_task
    wg.name = 'test_run_calcfunction'
    wg.run()
    print('state: ', wg.state)
    # print("results: ", results[])
    assert wg.tasks.sumdiff2.outputs.sum.value == 9


def test_dynamic_inputs() -> None:
    """Test dynamic inputs.
    For dynamic inputs, we allow the user to define the inputs manually.
    """

    @task.calcfunction()
    def add(**kwargs):
        return kwargs['x'] + kwargs['y']

    wg = WorkGraph('test_dynamic_inputs')
    wg.add_task(add, name='add1', x=orm.Int(1), y=orm.Int(2))
    # the top-level inputs are dynamic
    assert wg.tasks.add1.inputs._metadata.dynamic is True
    # the kwargs inputs are dynamic
    assert wg.tasks.add1.inputs.kwargs._metadata.dynamic is True
    assert wg.tasks.add1.inputs.kwargs._link_limit == 1e6
    wg.run()
    assert wg.tasks.add1.outputs.result.value == 3
    # the outputs are dynamic as well
    assert wg.tasks.add1.outputs._metadata.dynamic is True
