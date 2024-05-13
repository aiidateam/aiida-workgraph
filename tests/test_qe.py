import aiida
import numpy as np

aiida.load_profile()


def test_kpoint():
    """Run simple calcfunction."""
    from aiida_workgraph import WorkGraph

    wg = WorkGraph(name="test_kpoint")
    kpoint1 = wg.nodes.new("AiiDAKpoint", "kpoint1")
    kpoint1.set({"mesh": [2, 2, 2]})
    wg.submit(wait=True)
    assert kpoint1.outputs[0].value.get_kpoints_mesh() == ([2, 2, 2], [0, 0, 0])


def test_pw_parameters():
    """Run simple calcfunction."""
    from aiida_workgraph import WorkGraph

    wg = WorkGraph(name="test_pw_parameters")
    pw_parameters1 = wg.nodes.new("AiiDADict", "pw_parameters1")
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
    wg.submit(wait=True)
    assert pw_parameters1.node.get_dict() == paras


def test_structure(wg_structure_si):
    """Run simple calcfunction."""
    wg = wg_structure_si
    wg.name = "test_structure"
    print(wg.to_dict())
    wg.submit(wait=True)
    assert len(wg_structure_si.nodes["structure1"].node.get_ase()) == 2


def test_pw_pseudo(wg_structure_si):
    """Run simple calcfunction."""
    wg = wg_structure_si
    wg.name = "test_pw_pseudo"
    pw_pseudo1 = wg.nodes.new("AiiDAPWPseudo", "pw_pseudo1")
    wg.links.new(wg.nodes["structure1"].outputs[0], pw_pseudo1.inputs["structure"])
    wg.submit(wait=True)
    # assert pw_pseudo1.node


def test_pw_dos_projwfc(wg_structure_si):
    """Run simple calcfunction."""

    wg = wg_structure_si
    wg.name = "test_pw_dos_projwfc"
    #
    pw_relax1 = wg.nodes.new("AiiDAPW", "pw_relax1")
    metadata = {
        "options": {
            "stash": {},
            "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 2},
        }
    }
    pw_relax1.set({"metadata": metadata})
    #
    pw_code = wg.nodes.new("AiiDACode", "pw_code", label="qe-7.2-pw@localhost")
    #
    pw_parameters1 = wg.nodes.new("AiiDADict", "pw_parameters1")
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
    kpoint1 = wg.nodes.new("AiiDAKpoint", "kpoint1")
    kpoint1.set({"mesh": [2, 2, 2]})
    #
    pw_pseudo1 = wg.nodes.new("AiiDAPWPseudo", "pw_pseudo1")
    #
    #
    dos1 = wg.nodes.new("AiiDADos", "dos1")
    dos_code = wg.nodes.new("AiiDACode", "dos_code", label="qe-7.2-dos@localhost")
    dos_parameters1 = wg.nodes.new("AiiDADict", "dos_parameters1")
    dos1.set({"metadata": metadata})
    #
    projwfc1 = wg.nodes.new("AiiDAProjwfc", "projwfc1")
    projwfc_code = wg.nodes.new(
        "AiiDACode", "projwfc_code", label="qe-7.2-projwfc@localhost"
    )
    projwfc_parameters1 = wg.nodes.new("AiiDADict", "projwfc_parameters1")
    projwfc1.set({"metadata": metadata})
    #
    wg.links.new(wg.nodes["structure1"].outputs[0], pw_pseudo1.inputs["structure"])
    wg.links.new(wg.nodes["structure1"].outputs[0], pw_relax1.inputs["structure"])
    wg.links.new(pw_parameters1.outputs[0], pw_relax1.inputs["parameters"])
    wg.links.new(kpoint1.outputs[0], pw_relax1.inputs["kpoints"])
    wg.links.new(pw_pseudo1.outputs[0], pw_relax1.inputs["pseudos"])
    wg.links.new(pw_code.outputs[0], pw_relax1.inputs["code"])
    #
    wg.links.new(pw_relax1.outputs["remote_folder"], dos1.inputs["parent_folder"])
    wg.links.new(dos_code.outputs[0], dos1.inputs["code"])
    wg.links.new(dos_parameters1.outputs[0], dos1.inputs["parameters"])
    #
    wg.links.new(pw_relax1.outputs["remote_folder"], projwfc1.inputs["parent_folder"])
    wg.links.new(projwfc_code.outputs[0], projwfc1.inputs["code"])
    wg.links.new(projwfc_parameters1.outputs[0], projwfc1.inputs["parameters"])
    wg.submit(wait=True, timeout=150)
    assert wg.state == "FINISHED"
    assert np.isclose(
        wg.nodes["pw_relax1"].node.outputs.output_parameters["energy"], -305.9228430484
    )


def test_pw_relax_workchain(structure_si):
    """Run simple calcfunction."""
    from aiida_workgraph import worknode, WorkGraph
    from aiida.orm import Dict, KpointsData, load_code, load_group
    from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

    @worknode.calcfunction()
    def pw_parameters(paras, relax_type):
        paras1 = paras.clone()
        paras1["CONTROL"]["calculation"] = relax_type
        return paras1

    structure_si.store()
    #
    code = load_code("qe-7.2-pw@localhost")
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
    pseudo_family = load_group("SSSP/1.3/PBEsol/efficiency")
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

    wg = WorkGraph("test_pw_relax")
    # structure node
    wg.nodes.new("AiiDANode", "si", pk=structure_si.pk)
    # pw node
    pw_relax1 = wg.nodes.new(PwRelaxWorkChain, name="pw_relax1")
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
    paras_node = wg.nodes.new(
        pw_parameters, "parameters", paras=paras, relax_type="relax"
    )
    wg.links.new(wg.nodes["si"].outputs[0], pw_relax1.inputs["structure"])
    wg.links.new(paras_node.outputs[0], pw_relax1.inputs["base.pw.parameters"])
    wg.submit(wait=True, timeout=200)
    assert wg.state == "FINISHED"
    # print(wg.nodes["pw_relax1"].node.outputs.output_parameters["energy"])
    assert np.isclose(
        wg.nodes["pw_relax1"].node.outputs.output_parameters["energy"], -292.02237503211
    )
