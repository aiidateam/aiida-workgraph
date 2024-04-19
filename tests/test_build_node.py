from aiida_workgraph import build_node, Node


def test_calcjob():
    """Generate a node for test."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    ArithmeticAddNode = build_node(ArithmeticAddCalculation)
    assert issubclass(ArithmeticAddNode, Node)
    # build from path
    ArithmeticAddNode = build_node(
        "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"
    )
    assert issubclass(ArithmeticAddNode, Node)


def test_workchain():
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    MultiplyAddWorkNode = build_node(MultiplyAddWorkChain)
    assert issubclass(MultiplyAddWorkNode, Node)
    # build from path
    MultiplyAddWorkNode = build_node(
        "aiida.workflows.arithmetic.multiply_add.MultiplyAddWorkChain"
    )
    assert issubclass(MultiplyAddWorkNode, Node)


def test_calcfunction():
    """Generate a node for test."""
    from aiida.engine import calcfunction

    @calcfunction
    def add(x, y):
        """Calculate the square root of a number."""
        return x + y

    AddNode = build_node(add)
    assert issubclass(AddNode, Node)


def test_function():
    """Generate a node for test."""
    from scipy.linalg import norm

    AddNode = build_node(norm)
    assert issubclass(AddNode, Node)
