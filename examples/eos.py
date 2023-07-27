from aiida import load_profile
from aiida_worktree import build_node, WorkTree, node
from aiida.orm import Dict, KpointsData, StructureData, load_code, load_group
from ase.build import bulk
from aiida.engine import calcfunction

load_profile()

# register node
ndata = {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
pw_node = build_node(ndata)


@node()
@calcfunction
def scale_structure(structure, scale):
    atoms = structure.get_ase()
    atoms.set_cell(atoms.cell * scale, scale_atoms=True)
    return StructureData(ase=atoms)


# set link limit to a large value so that it can gather the result.
@node(
    inputs=[["General", "datas", {"link_limit": 100}]],
    outputs=[["General", "eos"]],
)
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


# ===============================================================================
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
# structure node
structure1 = nt.nodes.new("AiiDANode", "si", value=si)
eos1 = nt.nodes.new(eos, name="eos")
# pw node
scales = [0.95, 1.0, 1.05]
for i in range(len(scales)):
    pw1 = nt.nodes.new(pw_node, name=f"pw1_{i}")
    scale1 = nt.nodes.new(scale_structure, name=f"scale_{i}", scale=scales[i])
    pw1.set(
        {
            "code": code,
            "parameters": paras,
            "kpoints": kpoints,
            "pseudos": pseudos,
            "metadata": metadata,
        }
    )
    nt.links.new(structure1.outputs[0], scale1.inputs["structure"])
    nt.links.new(scale1.outputs[0], pw1.inputs["structure"])
    nt.links.new(pw1.outputs["output_parameters"], eos1.inputs[0])
nt.submit(wait=True, timeout=300)
print("eos: ", eos1.node.outputs.eos.get_dict())
