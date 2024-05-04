from aiida import load_profile
from aiida_workgraph import build_node, WorkGraph, node
from aiida.orm import Dict, KpointsData, StructureData, load_code, load_group
from ase.build import bulk

load_profile()

# register node
PwCalculation = build_node(
    {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
)


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


@node.graph_builder()
def eos_workgraph(structure=None, inputs=None, run_relax=True, scales=None):
    wg = WorkGraph()
    wg.context = {"current_structure": structure}
    # Load the pseudopotential family.
    pseudo_family = load_group("SSSP/1.3/PBEsol/efficiency")
    pseudos = pseudo_family.get_pseudos(structure=structure)
    # -------------relax----------------
    relax_node = wg.nodes.new(PwCalculation, name="relax1")
    relax_inputs = inputs.get("relax", {})
    relax_inputs["pseudos"] = pseudos
    relax_inputs["structure"] = "{{current_structure}}"
    relax_node.set(relax_inputs)
    # save the relaxed structure to the context
    relax_node.to_context = [["output_structure", "current_structure"]]
    # -------------eos euqatioin----------------
    eos1 = wg.nodes.new(eos, name="eos", datas="{{pw_result}}")
    eos1.wait = []
    # create pw node for each scaled structure
    for i in range(len(scales)):
        # -------------scale structure----------------
        # use the structure from the context
        scale1 = wg.nodes.new(
            scale_structure,
            name=f"scale_{i}",
            structure="{{structure}}",
            scale=scales[i],
        )
        scale1.wait = ["relax1"]
        # -------------scf----------------
        pw1 = wg.nodes.new(PwCalculation, name=f"pw_{i}")
        scf_inputs = inputs.get("scf", {})
        pw1.set(scf_inputs)
        # save the output_parameters to the context
        pw1.to_context = [["output_parameters", f"pw_result.s_{i}"]]
        # link the scale structure to the pw node
        wg.links.new(scale1.outputs[0], pw1.inputs["structure"])
        eos1.wait.append(pw1.name)
    if not run_relax:
        wg.nodes.delete("relax1")
    return wg


# ===============================================================================
# create input structure node
si = StructureData(ase=bulk("Si"))
# create the PW node
code = load_code("qe-7.2-pw@localhost")
kpoints = KpointsData()
kpoints.set_kpoints_mesh([2, 2, 2])
#
metadata = {
    "options": {
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        },
    }
}

inputs = {
    "relax": {
        "code": code,
        "parameters": Dict(
            {
                "CONTROL": {
                    "calculation": "vc-relax",
                },
                "SYSTEM": {
                    "ecutwfc": 30,
                    "ecutrho": 240,
                    "occupations": "smearing",
                    "smearing": "gaussian",
                    "degauss": 0.1,
                },
            }
        ),
        "kpoints": kpoints,
        "metadata": metadata,
    },
    "scf": {
        "code": code,
        "parameters": Dict(
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
        ),
        "kpoints": kpoints,
        "metadata": metadata,
    },
}

# -----------------------------------------------------------
wg = eos_workgraph(structure=si, inputs=inputs, scales=[0.95, 1.0, 1.05])
wg.submit(wait=True, timeout=300)
print("eos: ", wg.nodes["eos"].outputs["result"].value.get_dict())
