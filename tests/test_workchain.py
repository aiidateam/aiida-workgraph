import aiida
from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

aiida.load_profile()


def test_workchain(wg_workchain):
    """Submit simple calcjob."""
    wg = wg_workchain
    wg.name = "test_workchain"
    wg.submit(wait=True, timeout=100)
    # print("results: ", results[])
    assert wg.nodes["multiply_add2"].node.outputs.result == 17


def test_build_workchain_inputs_outputs():
    """Submit simple calcjob."""
    from aiida_workgraph import build_node

    node = build_node(MultiplyAddWorkChain)()
    assert len(node.inputs) == 10
    assert len(node.outputs) == 3


def test_build_workchain():
    """Submit simple calcjob."""
    from aiida.orm import load_code, Int
    from aiida_workgraph import WorkGraph

    code = load_code("add@localhost")
    wg = WorkGraph(name="test_debug_math")
    code1 = wg.nodes.new("AiiDACode", "code1", pk=code.pk)
    multiply_add1 = wg.nodes.new(
        MultiplyAddWorkChain,
        "multiply_add1",
        x=Int(4).store(),
        y=Int(2).store(),
        z=Int(3).store(),
    )
    multiply_add2 = wg.nodes.new(
        MultiplyAddWorkChain,
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    wg.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    wg.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    wg.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    wg.submit(wait=True, timeout=100)
    assert wg.nodes["multiply_add2"].node.outputs.result == 17
    # reload wg
    wg1 = WorkGraph.load(wg.pk)
    assert wg1.nodes["multiply_add2"].node.outputs.result == 17
