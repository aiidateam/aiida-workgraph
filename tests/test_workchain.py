import aiida

aiida.load_profile()


def test_workchain(wt_workchain):
    """Submit simple calcjob."""
    wt = wt_workchain
    wt.name = "test_workchain"
    wt.submit(wait=True, timeout=100)
    # print("results: ", results[])
    assert wt.nodes["multiply_add2"].node.outputs.result == 17


def test_build_workchain_inputs_outputs(build_workchain):
    """Submit simple calcjob."""

    node = build_workchain()
    assert len(node.inputs) == 10
    assert len(node.outputs) == 3


def test_build_workchain(build_workchain):
    """Submit simple calcjob."""
    from aiida.orm import load_code, Int
    from aiida_worktree import WorkTree

    code = load_code("add@localhost")
    wt = WorkTree(name="test_debug_math")
    code1 = wt.nodes.new("AiiDACode", "code1", value=code.pk)
    multiply_add1 = wt.nodes.new(
        build_workchain,
        "multiply_add1",
        x=Int(4).store(),
        y=Int(2).store(),
        z=Int(3).store(),
    )
    multiply_add2 = wt.nodes.new(
        build_workchain,
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    wt.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    wt.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    wt.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    wt.submit(wait=True, timeout=100)
    assert wt.nodes["multiply_add2"].node.outputs.result == 17
    # reload wt
    wt1 = WorkTree.load(wt.pk)
    assert wt1.nodes["multiply_add2"].node.outputs.result == 17
