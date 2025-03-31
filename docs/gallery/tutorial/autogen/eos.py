"""
==================================
Equation of state (EOS) WorkGraph
==================================

"""

# %%
# To run this tutorial, you need to install aiida-workgraph and set up a AiiDA profile.
#
# Create the calcfunction task
# ============================
#


from aiida import orm
from aiida_workgraph import task

#
# explicitly define the output socket name to match the return value of the function
@task.calcfunction(outputs=[{"name": "structures"}])
def scale_structure(structure, scales):
    """Scale the structure by the given scales."""
    atoms = structure.get_ase()
    structures = {}
    for i in range(len(scales)):
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scales[i], scale_atoms=True)
        structure = orm.StructureData(ase=atoms1)
        structures[f"s_{i}"] = structure
    return {"structures": structures}


@task.calcfunction()
# because this is a calcfunction, and the input datas are dynamic, we need use **datas.
def eos(**datas):
    """Fit the EOS of the data."""
    from ase.eos import EquationOfState

    #
    volumes = []
    energies = []
    for _, data in datas.items():
        volumes.append(data.dict.volume)
        energies.append(data.dict.energy)
        unit = data.dict.energy_units
    #
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    eos = orm.Dict({"unit": unit, "v0": v0, "e0": e0, "B": B})
    return eos


# %%
# Build the workgraph
# ===================
# Three steps:
#
# - create an empty WorkGraph
# - add tasks: scale_structure, scf and eos.
# - link the output and input sockets for the tasks.
#
#

from aiida_workgraph import WorkGraph, active_map_zone
from aiida_quantumespresso.calculations.pw import PwCalculation


def eos_workgraph(
    structure: orm.StructureData = None, scales: list = None, scf_inputs: dict = None
):
    with WorkGraph("eos_tutorial") as wg:
        wg.add_task(scale_structure, name="scale", structure=structure, scales=scales)
        with active_map_zone(wg.tasks.scale.outputs.structures) as map_zone:
            scf_task = map_zone.add_task(
                PwCalculation, name=f"scf", structure=map_zone.item
            )
            scf_task.set(scf_inputs)
        wg.add_task(eos, name="eos", datas=scf_task.outputs.output_parameters)
        return wg


# %%
# Prepare inputs and run
# ----------------------
# If you are running in a jupyter notebook, you can visualize the workgraph directly.
#


from aiida import load_profile
from aiida.common.exceptions import NotExistent
from aiida.orm import (
    Dict,
    load_code,
    load_group,
    InstalledCode,
    load_computer,
)
from ase.build import bulk

#
load_profile()
# create pw code
try:
    pw_code = load_code(
        "qe-7.2-pw@localhost"
    )  # The computer label can also be omitted here
except NotExistent:
    pw_code = InstalledCode(
        computer=load_computer("localhost"),
        filepath_executable="pw.x",
        label="qe-7.2-pw",
        default_calc_job_plugin="quantumespresso.pw",
    ).store()
#
si = orm.StructureData(ase=bulk("Si"))
pw_paras = Dict(
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
# Load the pseudopotential family.
pseudo_family = load_group("SSSP/1.3/PBEsol/efficiency")
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
#
kpoints = orm.KpointsData()
kpoints.set_kpoints_mesh([3, 3, 3])
pseudos = pseudo_family.get_pseudos(structure=si)
scf_inputs = {
    "code": pw_code,
    "parameters": pw_paras,
    "kpoints": kpoints,
    "pseudos": pseudos,
    "metadata": metadata,
}
# -------------------------------------------------------
# set the input parameters for each task
wg = eos_workgraph(si, [0.95, 1.0, 1.05], scf_inputs)
print("Waiting for the workgraph to finish...")
wg.submit(wait=True, timeout=300)
# one can also run the workgraph directly
# wg.run()
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg


# %%
# Print out the results:
#


data = wg.tasks.eos.outputs.result.value.get_dict()
print("B: {B}\nv0: {v0}\ne0: {e0}\nv0: {v0}".format(**data))

# %%
# Use it inside another workgraph
# ===============================
# In the above example, we call the `eos_workflow` directly to generate the workgraph and submit it. We can also use it as a task inside another workgraph. In this way, we can create a nested workflow.
#
# For example, we want to combine relax with eos.
#
# We can use the `graph_builder` decorator. The Graph Builder allow user to create a dynamic workflow based on the input value, as well as nested workflows.
#

from aiida_workgraph import WorkGraph, task, active_map_zone, active_graph


@task.graph_builder(outputs=[{"name": "result", "from": "eos.result"}])
def eos_workgraph(
    structure: orm.StructureData = None, scales: list = None, scf_inputs: dict = None
):
    with WorkGraph("eos") as wg:
        wg.add_task(scale_structure, name="scale", structure=structure, scales=scales)
        with active_map_zone(wg.tasks.scale.outputs.structures) as map_zone:
            scf_task = map_zone.add_task(
                PwCalculation, name=f"scf", structure=map_zone.item
            )
            scf_task.set(scf_inputs)
        wg.add_task(eos, name="eos", datas=scf_task.outputs.output_parameters)
        return wg


# %%
# Use it inside another workgraph
# -------------------------------
# For example, we want to combine relax with eos.
#


from aiida_workgraph import WorkGraph
from copy import deepcopy
from aiida_quantumespresso.calculations.pw import PwCalculation

# -------------------------------------------------------
relax_pw_paras = deepcopy(pw_paras)
relax_pw_paras["CONTROL"]["calculation"] = "vc-relax"
relax_inputs = {
    "structure": si,
    "code": pw_code,
    "parameters": relax_pw_paras,
    "kpoints": kpoints,
    "pseudos": pseudos,
    "metadata": metadata,
}
# -------------------------------------------------------
wg = WorkGraph("relax_eos")
wg.add_task(PwCalculation, name="relax")
wg.tasks.relax.set(relax_inputs)
eos_wg_task = wg.add_task(
    eos_workgraph,
    name="eos",
    structure=wg.tasks.relax.outputs.output_structure,
    scales=[0.95, 1.0, 1.05],
    scf_inputs=scf_inputs,
)
# -------------------------------------------------------
print("Waiting for the workgraph to finish...")
wg.run()
print(
    "\nResult: \nB: {B}\nv0: {v0}\ne0: {e0}\nv0: {v0}".format(
        **wg.tasks.eos.outputs.result.value.get_dict()
    )
)
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg

# %%
# Summary
# =======
# There are many ways to create the workflow using graph builder. For example, one can add the relax step inside the `eos_workgraph`, and add a `run_relax` argument to control the logic.
