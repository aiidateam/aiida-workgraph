from aiida_workgraph import build_node, WorkNode, WorkGraph


def test_calcjob():
    """Generate a node for test."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    # build from the class directly
    ArithmeticAddNode = build_node(ArithmeticAddCalculation)
    assert issubclass(ArithmeticAddNode, WorkNode)
    # build from path
    ArithmeticAddNode = build_node(
        "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"
    )
    assert issubclass(ArithmeticAddNode, WorkNode)
    # use the class directly
    wg = WorkGraph()
    add1 = wg.nodes.new(ArithmeticAddCalculation, name="add1")
    assert add1.name == "add1"


def test_workchain():
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    MultiplyAddWorkNode = build_node(MultiplyAddWorkChain)
    assert issubclass(MultiplyAddWorkNode, WorkNode)
    # build from path
    MultiplyAddWorkNode = build_node(
        "aiida.workflows.arithmetic.multiply_add.MultiplyAddWorkChain"
    )
    assert issubclass(MultiplyAddWorkNode, WorkNode)


def test_calcfunction():
    """Generate a node for test."""
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
    AddNode = build_node(add)
    assert issubclass(AddNode, WorkNode)
    # define outputs explicitly
    AddNode = build_node(
        add_minus, outputs=[["General", "sum"], ["General", "difference"]]
    )
    assert issubclass(AddNode, WorkNode)
    assert "sum" in AddNode().outputs.keys()
    # use the class directly
    wg = WorkGraph()
    add1 = wg.nodes.new(add, name="add1")
    assert "result" in add1.outputs.keys()
    assert add1.name == "add1"


def test_function():
    """Generate a node for test."""
    from scipy.linalg import norm

    AddNode = build_node(norm)
    assert issubclass(AddNode, WorkNode)
