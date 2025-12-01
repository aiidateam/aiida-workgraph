from aiida_workgraph import WorkGraph, task
from aiida_workgraph.task import TaskHandle


def test_calcjob():
    """Generate a task for test."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    # build from the class directly
    ArithmeticAddTask = task(ArithmeticAddCalculation)
    assert isinstance(ArithmeticAddTask, TaskHandle)
    # use the class directly
    wg = WorkGraph()
    add1 = wg.add_task(ArithmeticAddCalculation, name='add1')
    assert add1.name == 'add1'


def test_workchain():
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    MultiplyAddWorkTask = task(MultiplyAddWorkChain)
    assert isinstance(MultiplyAddWorkTask, TaskHandle)


def test_calcfunction():
    """Generate a task for test."""
    from aiida.engine import calcfunction

    @calcfunction
    def add(x, y):
        """Calculate the square root of a number."""
        return x + y

    @calcfunction
    def add_minus(x, y):
        """Calculate the square root of a number."""
        return {'sum': x + y, 'difference': x - y}

    # build from callable
    AddTask = task(add)
    assert isinstance(AddTask, TaskHandle)
    # define outputs explicitly
    AddTask = task(outputs=['sum', 'difference'])(add_minus)
    assert isinstance(AddTask, TaskHandle)
    assert 'sum' in AddTask()._task.get_output_names()
    # use the class directly
    wg = WorkGraph()
    add1 = wg.add_task(add, name='add1')
    assert 'result' in add1.get_output_names()
    assert add1.name == 'add1'

    AddTask_outputs_list = task(outputs=['sum', 'difference'])(add_minus)
    assert isinstance(AddTask_outputs_list, TaskHandle)
    assert 'sum' in AddTask_outputs_list()._task.get_output_names()


def test_function():
    """Generate a task for test."""
    from scipy.linalg import norm

    AddTask = task(norm)
    assert isinstance(AddTask, TaskHandle)
