from aiida_workgraph import build_task, Task, WorkGraph


def test_calcjob():
    """Generate a task for test."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    # build from the class directly
    ArithmeticAddTask = build_task(ArithmeticAddCalculation)
    assert issubclass(ArithmeticAddTask, Task)
    # build from path
    ArithmeticAddTask = build_task(
        "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"
    )
    assert issubclass(ArithmeticAddTask, Task)
    # use the class directly
    wg = WorkGraph()
    add1 = wg.add_task(ArithmeticAddCalculation, name="add1")
    assert add1.name == "add1"


def test_workchain():
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    MultiplyAddWorkTask = build_task(MultiplyAddWorkChain)
    assert issubclass(MultiplyAddWorkTask, Task)
    # build from path
    MultiplyAddWorkTask = build_task(
        "aiida.workflows.arithmetic.multiply_add.MultiplyAddWorkChain"
    )
    assert issubclass(MultiplyAddWorkTask, Task)


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
        return {"sum": x + y, "difference": x - y}

    # build from callable
    AddTask = build_task(add)
    assert issubclass(AddTask, Task)
    # define outputs explicitly
    AddTask = build_task(
        add_minus,
        outputs=[
            {"identifier": "Any", "name": "sum"},
            {"identifier": "Any", "name": "difference"},
        ],
    )
    assert issubclass(AddTask, Task)
    assert "sum" in AddTask().outputs.keys()
    # use the class directly
    wg = WorkGraph()
    add1 = wg.add_task(add, name="add1")
    assert "result" in add1.outputs.keys()
    assert add1.name == "add1"


def test_function():
    """Generate a task for test."""
    from scipy.linalg import norm

    AddTask = build_task(norm)
    assert issubclass(AddTask, Task)
