import aiida
import numpy as np

aiida.load_profile()


def test_pw_relax_workchain(structure_si):
    """Run simple calcfunction."""
    from aiida_worktree import build_node, WorkTree
    from aiida import orm

    # register node
    pw_relax_node = build_node(
        {"path": "aiida_quantumespresso.workflows.pw.relax.PwRelaxWorkChain"}
    )
    code = orm.load_code("pw-7.2@localhost")
    wt = WorkTree("test_pw_relax")
    pw_relax1 = wt.nodes.new(pw_relax_node, name="pw_relax1")
    pw_relax1.set_from_protocol(
        code, structure_si, protocol="fast", pseudo_family="SSSP/1.2/PBEsol/efficiency"
    )
    wt.submit(wait=True, timeout=200)
    assert wt.state == "FINISHED"
    # print(wt.nodes["pw_relax1"].node.outputs.output_parameters["energy"])
    assert np.isclose(
        wt.nodes["pw_relax1"].node.outputs.output_parameters["energy"], -308.46262827125
    )
