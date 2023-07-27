import aiida
import numpy as np

aiida.load_profile()


def test_kpoint():
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    nt = WorkTree(name="test_kpoint")
    kpoint1 = nt.nodes.new("AiiDAKpoint", "kpoint1")
    kpoint1.set({"mesh": [2, 2, 2]})
    nt.submit(wait=True)
    assert kpoint1.node.get_kpoints_mesh() == ([2, 2, 2], [0, 0, 0])


def test_pw_parameters():
    """Run simple calcfunction."""
    from aiida_worktree import WorkTree

    nt = WorkTree(name="test_pw_parameters")
    pw_parameters1 = nt.nodes.new("AiiDADict", "pw_parameters1")
    paras = {
        "CONTROL": {
            "calculation": "scf",
        },
        "SYSTEM": {
            "ecutwfc": 30,
            "ecutrho": 240,
        },
    }
    pw_parameters1.set({"value": paras})
    nt.submit(wait=True)
    assert pw_parameters1.node.get_dict() == paras


def test_structure(nt_structure_si):
    """Run simple calcfunction."""
    nt = nt_structure_si
    nt.name = "test_structure"
    nt.submit(wait=True)
    assert len(nt_structure_si.nodes["structure1"].node.get_ase()) == 2


def test_pw_pseudo(nt_structure_si):
    """Run simple calcfunction."""
    nt = nt_structure_si
    nt.name = "test_pw_pseudo"
    pw_pseudo1 = nt.nodes.new("AiiDAPWPseudo", "pw_pseudo1")
    nt.links.new(nt.nodes["structure1"].outputs[0], pw_pseudo1.inputs["structure"])
    nt.submit(wait=True)
    # assert pw_pseudo1.node


def test_pw_dos_projwfc(nt_structure_si):
    """Run simple calcfunction."""

    nt = nt_structure_si
    nt.name = "test_pw_dos_projwfc"
    #
    pw_relax1 = nt.nodes.new("AiiDAPW", "pw_relax1")
    metadata = {
        "options": {
            "stash": {},
            "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 2},
        }
    }
    pw_relax1.set({"metadata": metadata})
    #
    pw_code = nt.nodes.new("AiiDACode", "pw_code")
    pw_code.set({"value": "pw-7.2@localhost"})
    #
    pw_parameters1 = nt.nodes.new("AiiDADict", "pw_parameters1")
    paras = {
        "CONTROL": {
            "calculation": "scf",
        },
        "SYSTEM": {
            "ecutwfc": 30,
            "ecutrho": 240,
        },
    }
    pw_parameters1.set({"value": paras})
    #
    kpoint1 = nt.nodes.new("AiiDAKpoint", "kpoint1")
    kpoint1.set({"mesh": [2, 2, 2]})
    #
    pw_pseudo1 = nt.nodes.new("AiiDAPWPseudo", "pw_pseudo1")
    #
    #
    dos1 = nt.nodes.new("AiiDADos", "dos1")
    dos_code = nt.nodes.new("AiiDACode", "dos_code")
    dos_code.set({"value": "dos-7.2@localhost"})
    dos_parameters1 = nt.nodes.new("AiiDADict", "dos_parameters1")
    dos1.set({"metadata": metadata})
    #
    projwfc1 = nt.nodes.new("AiiDAProjwfc", "projwfc1")
    projwfc_code = nt.nodes.new("AiiDACode", "projwfc_code")
    projwfc_code.set({"value": "projwfc-7.2@localhost"})
    projwfc_parameters1 = nt.nodes.new("AiiDADict", "projwfc_parameters1")
    projwfc1.set({"metadata": metadata})
    #
    nt.links.new(nt.nodes["structure1"].outputs[0], pw_pseudo1.inputs["structure"])
    nt.links.new(nt.nodes["structure1"].outputs[0], pw_relax1.inputs["structure"])
    nt.links.new(pw_parameters1.outputs[0], pw_relax1.inputs["parameters"])
    nt.links.new(kpoint1.outputs[0], pw_relax1.inputs["kpoints"])
    nt.links.new(pw_pseudo1.outputs[0], pw_relax1.inputs["pseudos"])
    nt.links.new(pw_code.outputs[0], pw_relax1.inputs["code"])
    #
    nt.links.new(pw_relax1.outputs["remote_folder"], dos1.inputs["parent_folder"])
    nt.links.new(dos_code.outputs[0], dos1.inputs["code"])
    nt.links.new(dos_parameters1.outputs[0], dos1.inputs["parameters"])
    #
    nt.links.new(pw_relax1.outputs["remote_folder"], projwfc1.inputs["parent_folder"])
    nt.links.new(projwfc_code.outputs[0], projwfc1.inputs["code"])
    nt.links.new(projwfc_parameters1.outputs[0], projwfc1.inputs["parameters"])
    nt.submit(wait=True, timeout=150)
    assert nt.state == "FINISHED"
    assert np.isclose(
        nt.nodes["pw_relax1"].node.outputs.output_parameters["energy"], -305.9228430484
    )


def test_pw_relax_workchain(structure_si):
    """Run simple calcfunction."""
    from aiida_worktree import build_node, node, WorkTree
    from aiida.engine import calcfunction
    from aiida.orm import Dict, KpointsData, load_code, load_group

    # register node
    ndata = {"path": "aiida_quantumespresso.workflows.pw.relax.PwRelaxWorkChain"}
    pw_relax_node = build_node(ndata)
    #
    @node()
    @calcfunction
    def pw_parameters(paras, relax_type):
        paras1 = paras.clone()
        paras1["CONTROL"]["calculation"] = relax_type
        return paras1

    #
    code = load_code("pw-7.2@localhost")
    paras = Dict(
        {
            "CONTROL": {
                "calculation": "scf",
            },
            "SYSTEM": {
                "ecutwfc": 30,
                "ecutrho": 240,
                "occupations": "smearing",
                "smearing": "gaussian",
                "degauss": 0.1,
            },
        }
    )
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([1, 1, 1])
    # Load the pseudopotential family.
    pseudo_family = load_group("SSSP/1.2/PBEsol/efficiency")
    pseudos = pseudo_family.get_pseudos(structure=structure_si)
    #
    metadata = {
        "options": {
            "resources": {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            },
        }
    }

    nt = WorkTree("test_pw_relax")
    # structure node
    nt.nodes.new("AiiDANode", "si", value=structure_si)
    # pw node
    pw_relax1 = nt.nodes.new(pw_relax_node, name="pw_relax1")
    pw_relax1.set(
        {
            "base": {
                "pw": {
                    "code": code,
                    "pseudos": pseudos,
                    "metadata": metadata,
                },
                "kpoints": kpoints,
            },
        }
    )
    paras_node = nt.nodes.new(
        pw_parameters, "parameters", paras=paras, relax_type="relax"
    )
    nt.links.new(nt.nodes["si"].outputs[0], pw_relax1.inputs["structure"])
    nt.links.new(paras_node.outputs[0], pw_relax1.inputs["base.pw.parameters"])
    nt.submit(wait=True, timeout=200)
    assert nt.state == "FINISHED"
    # print(nt.nodes["pw_relax1"].node.outputs.output_parameters["energy"])
    assert np.isclose(
        nt.nodes["pw_relax1"].node.outputs.output_parameters["energy"], -292.02237503211
    )
