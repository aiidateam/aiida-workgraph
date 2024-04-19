from aiida_workgraph import node, WorkGraph, build_node
from aiida.orm import Int, load_code
from aiida import load_profile

load_profile()


arithmetic_add = build_node(
    {"path": "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"}
)


@node.calcfunction()
def decorated_add(x, y, t=1):
    import time

    time.sleep(t.value)
    return x + y


@node.calcfunction()
def decorated_multiply(x, y, t=1):
    import time

    time.sleep(t.value)
    return x * y


@node.group(outputs=[["multiply.result", "result"]])
def add_multiply_group(x, y, z, t=2):
    wt = WorkGraph("add_multiply_group")
    add1 = wt.nodes.new(decorated_add, name="add1", x=x, y=y, t=t)
    multiply = wt.nodes.new(decorated_multiply, name="multiply", x=z, t=t)
    # link the output of int node to the input of add node
    wt.links.new(add1.outputs[0], multiply.inputs["y"])
    return wt


"""Use to test the engine."""
code = load_code("add@localhost")
x = Int(2)
wt = WorkGraph(name="test_run_order")
adds = []
add0 = wt.nodes.new(decorated_add, "add_0", x=x, y=x, t=Int(1))
multiply0 = wt.nodes.new(decorated_multiply, "multiply_0", x=x, y=x, t=Int(1))
add_multiply_group0 = wt.nodes.new(
    add_multiply_group, "add_multiply_group_0", x=x, y=x, z=x, t=Int(1)
)
arithmetic_add0 = wt.nodes.new(arithmetic_add, "arithmetic_add_0", x=x, y=x, code=code)
arithmetic_add1 = wt.nodes.new(arithmetic_add, "arithmetic_add_1", x=x, y=x, code=code)
arithmetic_add0.set({"metadata.options.sleep": 15})
arithmetic_add1.set({"metadata.options.sleep": 1})
add1 = wt.nodes.new(decorated_multiply, "add_1", x=x, y=x)
wt.links.new(add0.outputs[0], multiply0.inputs["x"])
wt.links.new(multiply0.outputs[0], add_multiply_group0.inputs["x"])
wt.links.new(arithmetic_add0.outputs["sum"], arithmetic_add1.inputs["x"])
wt.links.new(add_multiply_group0.outputs[0], add1.inputs["x"])
wt.links.new(arithmetic_add1.outputs["sum"], add1.inputs["y"])
wt.submit()
