from aiida_workgraph import WorkGraph, Task, task
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


def test_error_handlers(add_code):
    """Test error handlers."""

    def handle_negative_sum(task: Task):
        """Handle negative sum by resetting the task and changing the inputs.
        self is the WorkGraph instance, thus we can access the tasks and the context.
        """
        # modify task inputs
        task.set_inputs(
            {
                'x': orm.Int(abs(task.inputs.x.value)),
                'y': orm.Int(abs(task.inputs['y'].value)),
            }
        )
        msg = 'Run error handler: handle_negative_sum.'
        return msg

    wg = WorkGraph('restart_graph')
    add1 = wg.add_task(ArithmeticAddCalculation, name='add1')
    add1.add_error_handler(
        {
            'handle_negative_sum': {
                'executor': handle_negative_sum,
                'exit_codes': [
                    ArithmeticAddCalculation.exit_codes.ERROR_NEGATIVE_NUMBER.status,
                ],
                'max_retries': 5,
                'kwargs': {},
            }
        }
    )
    assert len(wg.tasks.add1.error_handlers) == 1
    wg.run(
        inputs={
            'add1': {'code': add_code, 'x': orm.Int(1), 'y': orm.Int(-2)},
        },
    )
    assert wg.tasks.add1.outputs.sum.value == 3


def test_error_handlers_graph_inputs(add_code):
    """Test error handlers with graph inputs linked to task inputs."""
    from aiida_workgraph.tasks.tests import BaseAddTask

    @task.graph()
    def restart_graph(code, x, y):
        BaseAddTask(code=code, x=x, y=y)

    g = restart_graph.build(code=add_code, x=1, y=-2)
    assert len(g.tasks.ArithmeticAddCalculation.error_handlers) == 1
    g1 = WorkGraph.from_dict(g.to_dict())
    assert len(g1.tasks.ArithmeticAddCalculation.error_handlers) == 1
    g.run()
    assert g.tasks.ArithmeticAddCalculation.outputs.sum.value == 3
