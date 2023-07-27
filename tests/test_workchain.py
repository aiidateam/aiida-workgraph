import aiida

aiida.load_profile()


def test_workchain(nt_workchain):
    """Submit simple calcjob."""
    nt = nt_workchain
    nt.name = "test_workchain"
    nt.submit(wait=True, timeout=100)
    # print("results: ", results[])
    assert nt.nodes["multiply_add2"].node.outputs.result == 17


def test_build_workchain_inputs_outputs(build_workchain):
    """Submit simple calcjob."""

    node = build_workchain()
    assert len(node.inputs) == 9
    assert len(node.outputs) == 1


def test_build_workchain(build_workchain):
    """Submit simple calcjob."""
    from aiida.orm import load_code, Int
    from aiida_worktree import WorkTree

    code = load_code("add@localhost")
    nt = WorkTree(name="test_debug_math")
    code1 = nt.nodes.new("AiiDACode", "code1", value=code.pk)
    multiply_add1 = nt.nodes.new(
        build_workchain,
        "multiply_add1",
        x=Int(4).store(),
        y=Int(2).store(),
        z=Int(3).store(),
    )
    multiply_add2 = nt.nodes.new(
        build_workchain,
        "multiply_add2",
        x=Int(2).store(),
        y=Int(3).store(),
    )
    nt.links.new(code1.outputs[0], multiply_add1.inputs["code"])
    nt.links.new(code1.outputs[0], multiply_add2.inputs["code"])
    nt.links.new(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    nt.submit(wait=True, timeout=100)
    assert nt.nodes["multiply_add2"].node.outputs.result == 17
