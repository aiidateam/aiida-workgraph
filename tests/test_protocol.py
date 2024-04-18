import aiida
import numpy as np

aiida.load_profile()


def test_pw_relax_protocol(structure_si):
    """Run simple calcfunction."""
    from aiida_worktree import build_node, WorkTree
    from aiida import orm

    # register node
    pw_relax_node = build_node(
        "aiida_quantumespresso.workflows.pw.relax.PwRelaxWorkChain"
    )
    code = orm.load_code("qe-7.2-pw@localhost")
    wt = WorkTree("test_pw_relax")
    pw_relax1 = wt.nodes.new(pw_relax_node, name="pw_relax1")
    pw_relax1.set_from_protocol(
        code,
        structure_si,
        protocol="fast",
    )
    wt.submit(wait=True, timeout=200)
    assert wt.state == "FINISHED"
    assert "base_final_scf" in wt.nodes["pw_relax1"].node.inputs
    energy = wt.nodes["pw_relax1"].node.outputs.output_parameters["energy"]
    assert np.isclose(energy, -308.46, atol=0.1)


def test_pw_relax_protocol_pop(structure_si):
    """Run simple calcfunction."""
    from aiida_worktree import build_node, WorkTree
    from aiida import orm

    # register node
    pw_relax_node = build_node(
        "aiida_quantumespresso.workflows.pw.relax.PwRelaxWorkChain"
    )
    code = orm.load_code("qe-7.2-pw@localhost")
    wt = WorkTree("test_pw_relax")
    pw_relax1 = wt.nodes.new(pw_relax_node, name="pw_relax1")
    pw_relax1.set_from_protocol(
        code,
        structure_si,
        protocol="fast",
    )
    # do not run "base_final_scf"
    pw_relax1.inputs["base_final_scf"].value = None
    wt.submit(wait=True, timeout=200)
    assert "base_final_scf" not in wt.nodes["pw_relax1"].node.inputs
    energy = wt.nodes["pw_relax1"].node.outputs.output_parameters["energy"]
    assert np.isclose(energy, -308.46, atol=0.1)
