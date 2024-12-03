"""
==================================
Equation of state (EOS) WorkGraph
==================================

"""

# %%
# To run this tutorial, you need to install aiida-workgraph and restart the daemon. Open a terminal and run:
#
# .. code-block:: console
#
#    pip install aiida-workgraph aiida-quantumespresso
#
# Restart (or start) the AiiDA daemon if needed:
#
# .. code-block:: console
#
#    verdi daemon restart
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


#
# Output result from context to the output socket
@task.graph_builder(outputs=[{"name": "result", "from": "context.result"}])
def all_scf(structures, scf_inputs):
    """Run the scf calculation for each structure."""
    from aiida_workgraph import WorkGraph
    from aiida_quantumespresso.calculations.pw import PwCalculation

    wg = WorkGraph()
    for key, structure in structures.items():
        pw1 = wg.add_task(PwCalculation, name=f"pw1_{key}", structure=structure)
        pw1.set(scf_inputs)
        # save the output parameters to the context
        pw1.set_context({f"result.{key}": "output_parameters"})
    return wg


#


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
# - add tasks: scale_structure, all_scf and eos.
# - link the output and input sockets for the tasks.
#
# Visualize the workgraph
# -----------------------
# If you are running in a jupyter notebook, you can visualize the workgraph directly.
#

from aiida_workgraph import WorkGraph

#
wg = WorkGraph("eos")
scale_structure1 = wg.add_task(scale_structure, name="scale_structure1")
all_scf1 = wg.add_task(all_scf, name="all_scf1")
eos1 = wg.add_task(eos, name="eos1")
wg.add_link(scale_structure1.outputs["structures"], all_scf1.inputs["structures"])
wg.add_link(all_scf1.outputs["result"], eos1.inputs["datas"])
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg


# %%
# Prepare inputs and run
# ----------------------
#


from aiida import load_profile
from aiida.common.exceptions import NotExistent
from aiida.orm import (
    Dict,
    KpointsData,
    StructureData,
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
wg.tasks["scale_structure1"].set({"structure": si, "scales": [0.95, 1.0, 1.05]})
wg.tasks["all_scf1"].set({"scf_inputs": scf_inputs})
print("Waiting for the workgraph to finish...")
wg.submit(wait=True, timeout=300)
# one can also run the workgraph directly
# wg.run()


# %%
# Print out the results:
#


data = wg.tasks["eos1"].outputs["result"].value.get_dict()
print("B: {B}\nv0: {v0}\ne0: {e0}\nv0: {v0}".format(**data))

# %%
# Use graph builder
# =================
# The Graph Builder allow user to create a dynamic workflow based on the input value, as well as nested workflows.
#

from aiida_workgraph import WorkGraph, task

#
@task.graph_builder(outputs=[{"name": "result", "from": "eos1.result"}])
def eos_workgraph(structure=None, scales=None, scf_inputs=None):
    wg = WorkGraph("eos")
    scale_structure1 = wg.add_task(
        scale_structure, name="scale_structure1", structure=structure, scales=scales
    )
    all_scf1 = wg.add_task(all_scf, name="all_scf1", scf_inputs=scf_inputs)
    eos1 = wg.add_task(eos, name="eos1")
    wg.add_link(scale_structure1.outputs["structures"], all_scf1.inputs["structures"])
    wg.add_link(all_scf1.outputs["result"], eos1.inputs["datas"])
    return wg


# %%
# Then we can use the `eos_workgraph` in two ways:
#
# - Direct run the function and generate the workgraph, then submit
# - Use it as a task inside another workgraph to create nested workflow.
#
# Use the graph builder directly
# ------------------------------
#

wg = eos_workgraph(structure=si, scales=[0.95, 1.0, 1.05], scf_inputs=scf_inputs)
# One can submit the workgraph directly
# wg.submit(wait=True, timeout=300)
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg

# %%
# Use it inside another workgraph
# -------------------------------
# For example, we want to combine relax with eos.
#


from aiida_workgraph import WorkGraph
from copy import deepcopy
from aiida_quantumespresso.calculations.pw import PwCalculation

#
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
relax_task = wg.add_task(PwCalculation, name="relax1")
relax_task.set(relax_inputs)
eos_wg_task = wg.add_task(
    eos_workgraph, name="eos1", scales=[0.95, 1.0, 1.05], scf_inputs=scf_inputs
)
wg.add_link(relax_task.outputs["output_structure"], eos_wg_task.inputs["structure"])
# -------------------------------------------------------
# One can submit the workgraph directly
# wg.submit(wait=True, timeout=300)

wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg

# %%
# Summary
# =======
# There are many ways to create the workflow using graph builder. For example, one can add the relax step inside the `eos_workgraph`, and add a `run_relax` argument to control the logic.
