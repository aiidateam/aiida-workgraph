from aiida import load_profile
from aiida_workgraph import build_node, WorkGraph, node
from aiida.orm import Dict, KpointsData, StructureData, load_code, load_group
from ase.build import bulk

load_profile()

# register node
ndata = {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
pw_node = build_node(ndata)


@node.calcfunction()
def scale_structure(structure, scale):
    atoms = structure.get_ase()
    atoms.set_cell(atoms.cell * scale, scale_atoms=True)
    return StructureData(ase=atoms)


@node.calcfunction()
# because this is a calcfunction, and the input datas are dynamic, we need use **datas.
def eos(**datas):
    from ase.eos import EquationOfState

    volumes = []
    energies = []
    for _key, data in datas.items():
        # it only gather the uuid of the data, so we need to load it.
        volumes.append(data.dict.volume)
        energies.append(data.dict.energy)
        unit = data.dict.energy_units
    #
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    eos = Dict(
        {
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

wt = WorkGraph("eos")
# structure node
structure1 = wt.nodes.new("AiiDANode", "si", value=si)
# get the result of each pw node from the context
eos1 = wt.nodes.new(eos, name="eos", datas="{{pw_result}}")
# create pw node for each scale
scales = [0.95, 1.0, 1.05]
for i in range(len(scales)):
    pw1 = wt.nodes.new(pw_node, name=f"pw1_{i}")
    scale1 = wt.nodes.new(scale_structure, name=f"scale_{i}", scale=scales[i])
    pw1.set(
        {
            "code": code,
            "parameters": paras,
            "kpoints": kpoints,
            "pseudos": pseudos,
            "metadata": metadata,
        }
    )
    pw1.to_ctx = [["output_parameters", f"pw_result.s_{i}"]]
    wt.links.new(structure1.outputs[0], scale1.inputs["structure"])
    wt.links.new(scale1.outputs[0], pw1.inputs["structure"])
    wt.ctrl_links.new(pw1.ctrl_outputs[0], eos1.ctrl_inputs[0])
wt.submit(wait=True, timeout=300)
print("eos: ", eos1.node.outputs.eos.get_dict())
