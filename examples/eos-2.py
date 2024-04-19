from aiida import load_profile, orm
from aiida_workgraph import build_node, WorkGraph, node
from ase.build import bulk

load_profile()


@node.calcfunction(outputs=[["General", "structures"]])
def scale_structure(structure, scales):
    atoms = structure.get_ase()
    structures = {}
    for i in range(len(scales)):
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scales[i], scale_atoms=True)
        structure = orm.StructureData(ase=atoms1)
        structures[f"s_{i}"] = structure
    return {"structures": structures}


# Output result from context
@node.group(outputs=[["ctx.result", "result"]])
def all_scf(structures, code, parameters, kpoints, pseudos, metadata):
    from aiida_workgraph import WorkGraph

    # register node
    ndata = {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
    pw_node = build_node(ndata)
    wt = WorkGraph("all_scf")
    # pw node
    for key, structure in structures.items():
        pw1 = wt.nodes.new(pw_node, name=f"pw1_{key}")
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
        pw1.to_ctx = [["output_parameters", f"result.{key}"]]
    return wt


@node.calcfunction()
# because this is a calcfunction, and the input datas are dynamic, we need use **datas.
def eos(**datas):
    from ase.eos import EquationOfState

    volumes = []
    energies = []
    for _key, data in datas.items():
        volumes.append(data.dict.volume)
        energies.append(data.dict.energy)
        unit = data.dict.energy_units
    #
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    eos = orm.Dict(
        {
            "unit": unit,
            "v0": v0,
            "e0": e0,
            "B": B,
        }
    )
    return eos


# ===================================================
# create input structure node
si = orm.StructureData(ase=bulk("Si"))
# create the PW node
code = orm.load_code("qe-7.2-pw@localhost")
paras = orm.Dict(
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
kpoints = orm.KpointsData()
kpoints.set_kpoints_mesh([2, 2, 2])
# Load the pseudopotential family.
pseudo_family = orm.load_group("SSSP/1.2/PBEsol/efficiency")
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
scale_structure1 = wt.nodes.new(
    scale_structure, name="scale_structure1", structure=si, scales=[0.95, 1.0, 1.05]
)
all_scf1 = wt.nodes.new(all_scf, name="all_scf1")
all_scf1.set(
    {
        "code": code,
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "metadata": metadata,
    }
)
eos1 = wt.nodes.new(eos, name="eos1")
wt.links.new(scale_structure1.outputs["structures"], all_scf1.inputs["structures"])
wt.links.new(all_scf1.outputs["result"], eos1.inputs["datas"])
wt.submit(wait=True, timeout=300)
