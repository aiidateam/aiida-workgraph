"""
================================
Computational materials science
================================

"""

# %%
# Introduction
# ============
# In this tutorial, you will use `AiiDA-WorkGraph` to carry out a DFT calculation using Quantum ESPRESSO.
#
# Requirements
# ------------
# To run this tutorial, you need to install `aiida-workgraph`, `aiida-quantumespresso` and `aiida-pseudo`. Open a terminal and run:
#
# .. code-block:: console
#
#    pip install aiida-workgraph aiida-quantumespresso aiida-pseudo
#    aiida-pseudo install sssp -x PBEsol
#
# Start the AiiDA daemon if needed:
#
# .. code-block:: console
#
#    verdi daemon start
#
# Start the web server
# --------------------
#
# Open a terminal, and run:
#
# .. code-block:: console
#
#    workgraph web start
#
# Then visit the page `http://127.0.0.1:8000/workgraph`, where you can view the workgraph later.
#
# Load the AiiDA profile.
#


from aiida import load_profile

load_profile()
#
# %%
# First workflow: calculate the energy of N2 molecule
# ===================================================
# Define a workgraph
# -------------------
# aiida-quantumespresso provides a CalcJob: `PwCalculation` to run a PW calculation. we can use it directly in the WorkGraph. The inputs and outputs of the task is automatically generated based on the `PwCalculation` CalcJob.
#

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_workgraph import WorkGraph

#
wg = WorkGraph("energy_n2")
pw1 = wg.add_task(PwCalculation, name="pw1")
pw1.to_html()
#
# visualize the task in jupyter-notebook
# pw1
#

# %%
# Prepare the inputs and submit the workflow
# ------------------------------------------
#
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
from ase.build import molecule

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
# create input structure
mol = molecule("N2")
mol.center(vacuum=1.5)
mol.pbc = True
structure_n2 = StructureData(ase=mol)
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
pseudos = pseudo_family.get_pseudos(structure=structure_n2)
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
# ------------------------- Set the inputs -------------------------
pw1.set(
    {
        "code": pw_code,
        "structure": structure_n2,
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "metadata": metadata,
    }
)
# ------------------------- Submit the calculation -------------------------
wg.submit(wait=True, timeout=200)
# ------------------------- Print the output -------------------------
print(
    "Energy of an un-relaxed N2 molecule: {:0.3f}".format(
        pw1.outputs["output_parameters"].value.get_dict()["energy"]
    )
)
#

# %%
# Generate node graph from the AiiDA process:
#

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Second workflow: atomization energy of N2 molecule
# ==================================================
#
# The atomization energy of :math:`N_2` is defined as the energy difference between the :math:`N_2` molecule and two isolated N atoms.
#
# .. code-block:: python
#
#    e_atomization = 2 * e_atom - e_molecule

# Define a calcfunction to calculate the atomization energy
# ---------------------------------------------------------
#

from aiida_workgraph import task

#
@task.calcfunction()
def atomization_energy(output_atom, output_mol):
    from aiida.orm import Float

    e = output_atom["energy"] * output_mol["number_of_atoms"] - output_mol["energy"]
    return Float(e)


# %%
# Create the structure of nitrogen Atom.
#

from ase import Atoms
from aiida.orm import StructureData

#
atoms = Atoms("N")
atoms.center(vacuum=1.5)
atoms.pbc = True
structure_n = StructureData(ase=atoms)

# %%
# Create a workgraph
# ------------------


from aiida_workgraph import WorkGraph
from aiida.orm import load_code

#
# load the PW code
pw_code = load_code("qe-7.2-pw@localhost")
#
wg = WorkGraph("atomization_energy")
#
# create the PW task
pw_n = wg.add_task(PwCalculation, name="pw_n")
pw_n.set(
    {
        "code": pw_code,
        "structure": structure_n,
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "metadata": metadata,
    }
)
pw_n2 = wg.add_task(PwCalculation, name="pw_n2")
pw_n2.set(
    {
        "code": pw_code,
        "structure": structure_n2,
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "metadata": metadata,
    }
)
# create the task to calculate the atomization energy
atomization = wg.add_task(atomization_energy, name="atomization_energy")
wg.add_link(pw_n.outputs["output_parameters"], atomization.inputs["output_atom"])
wg.add_link(pw_n2.outputs["output_parameters"], atomization.inputs["output_mol"])
wg.to_html()


# %%
# Submit the workgraph and print the atomization energy.
#


wg.submit(wait=True, timeout=300)
print(
    "Atomization energy: {:0.3f} eV".format(atomization.outputs["result"].value.value)
)


# %%
# If you start the web app (`workgraph web start`), you can visit the page http://127.0.0.1:8000/workgraph to view the tasks.
#
# You can also generate node graph from the AiiDA process:
#


from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Use already existing workchain
# ===============================
# Can we register a task from a workchain? Can we set the a input item of a namespace? Yes, we can!
#
# In the `PwRelaxWorkChain`, one can set the relax type (`calculation` key) in the input namespace `base.pw.parameters`. Now we create a new task to update the pw parameters.
#

from aiida_workgraph import task


@task.calcfunction()
def pw_parameters(paras, relax_type):
    paras1 = paras.clone()
    paras1["CONTROL"]["calculation"] = relax_type
    return paras1


# %%
# Now, we create the workgraph to relax the structure of N2 molecule.
#

from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

#
wg = WorkGraph("test_pw_relax")
# pw task
pw_relax1 = wg.add_task(PwRelaxWorkChain, name="pw_relax1")
# Load the pseudopotential family.
pseudos = pseudo_family.get_pseudos(structure=structure_n2)
pw_relax1.set(
    {
        "base": {
            "pw": {"code": pw_code, "pseudos": pseudos, "metadata": metadata},
            "kpoints": kpoints,
        },
        "structure": structure_n2,
    },
)
paras_task = wg.add_task(pw_parameters, "parameters", paras=paras, relax_type="relax")
wg.add_link(paras_task.outputs[0], pw_relax1.inputs["base.pw.parameters"])
# One can submit the workgraph directly
# wg.submit(wait=True, timeout=200)
# print(
#     "\nEnergy of a relaxed N2 molecule: {:0.3f}".format(
#         pw_relax1.node.outputs.output_parameters.get_dict()["energy"]
#     )
# )


# %%
# Use `protocol` to set input parameters (Experimental)
# ====================================================
# The aiida-quantumespresso package supports setting input parameters from protocol. For example, the PwRelaxWorkChain has a `get_builder_from_protocol` method. In this tutorial, we will show how to use the `protocol` to set the input parameters inside the WorkGraph.
#

from aiida_workgraph import build_task, WorkGraph
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain
from ase.build import bulk
from aiida import orm
from pprint import pprint

#
pw_code = orm.load_code("qe-7.2-pw@localhost")
wg = WorkGraph("test_pw_relax")
structure_si = orm.StructureData(ase=bulk("Si"))
pw_relax1 = wg.add_task(PwRelaxWorkChain, name="pw_relax1")
# set the inputs from the protocol
# this will call the `PwRelaxWorkChain.get_builder_from_protocol` method
# to set the inputs of the workchain
pw_relax1.set_from_protocol(
    pw_code, structure_si, protocol="fast", pseudo_family="SSSP/1.2/PBEsol/efficiency"
)
# we can now inspect the inputs of the workchain
print("The inputs for the PwBaseWorkchain are:")
print("-" * 80)
pprint(pw_relax1.inputs["base"].value)
print("\nThe input parameters for pw are:")
print("-" * 80)
pprint(pw_relax1.inputs["base"].value["pw"]["parameters"].get_dict())


# %%
# One can also adjust the parameters of the `PwRelaxWorkChain` to from protocol.
#

# For example, we want to remove the `base_final_scf` from the inputs, so that the `PwRelaxWorkChain` will not run the `base_final_scf` step.
pw_relax1.inputs["base_final_scf"].value = None
# submit the workgraph
# wg.submit(wait=True, timeout=200)
