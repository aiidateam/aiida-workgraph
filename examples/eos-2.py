from aiida import load_profile
from aiida_worktree import build_node, WorkTree, node
from aiida.orm import (
    Dict,
    KpointsData,
    StructureData,
    load_code,
    load_group,
    List,
    load_node,
)
from ase.build import bulk
from aiida.engine import calcfunction

load_profile()


@node()
@calcfunction
def scale_structure(structure, scales):

    atoms = structure.get_ase()
    structures = []
    for scale in scales:
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scale, scale_atoms=True)
        structure = StructureData(ase=atoms1).store()
        structures.append(structure.uuid)
    return List(list=structures)


# the structure_uuids is used to generate the worktree dynamically.
@node.group(outputs=[["gather1", "result", "result"]])
def all_scf(structure_uuids, code, parameters, kpoints, pseudos, metadata):
    # register node
    ndata = {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
    pw_node = build_node(ndata)

    nt = WorkTree("all_scf")
    gather1 = nt.nodes.new("AiiDAGather", name="gather1")
    # pw node
    for i in range(len(structure_uuids)):
        structure = load_node(structure_uuids[i])
        pw1 = nt.nodes.new(pw_node, name=f"pw1_{i}")
        pw1.set(
            {
                "code": code,
                "parameters": parameters,
                "kpoints": kpoints,
                "pseudos": pseudos,
                "metadata": metadata,
                "structure": structure,
            }
        )
        nt.links.new(pw1.outputs["output_parameters"], gather1.inputs[0])
    return nt


# set link limit to a large value so that it can gather the result.
@node()
@calcfunction
def eos(datas):
    from ase.eos import EquationOfState
    from aiida.orm import load_node

    volumes = []
    energies = []
    for data in datas:
        # it only gather the uuid of the data, so we need to load it.
        data = load_node(data)
        volumes.append(data.dict.volume)
        energies.append(data.dict.energy)
        unit = data.dict.energy_units
    #
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    eos = Dict(
        {
            "volumes": volumes,
            "energies": energies,
            "unit": unit,
            "v0": v0,
            "e0": e0,
            "B": B,
        }
    )
    return eos


# ===================================================
# create input structure node
si = StructureData(ase=bulk("Si"))
# create the PW node
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
kpoints.set_kpoints_mesh([2, 2, 2])
# Load the pseudopotential family.
pseudo_family = load_group("SSSP/1.2/PBEsol/efficiency")
pseudos = pseudo_family.get_pseudos(structure=si)
#
metadata = {
    "options": {
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        },
    }
}

nt = WorkTree("eos")
scale_structure1 = nt.nodes.new(
    scale_structure, name="scale_structure1", structure=si, scales=[0.95, 1.0, 1.05]
)
all_scf1 = nt.nodes.new(all_scf, name="all_scf1")
all_scf1.set(
    {
        "code": code,
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "metadata": metadata,
    }
)
eos1 = nt.nodes.new(eos, name="eos1")
nt.links.new(scale_structure1.outputs["result"], all_scf1.inputs["structure_uuids"])
nt.links.new(all_scf1.outputs["result"], eos1.inputs["datas"])

nt.submit(wait=True, timeout=300)
