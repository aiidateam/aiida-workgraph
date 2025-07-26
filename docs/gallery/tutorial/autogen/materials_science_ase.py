"""
================================
Computational materials science
================================

"""

# %%
# Introduction
# ============
# In this tutorial, you will use `AiiDA-WorkGraph` to carry out a material science calculation with using ASE package.
#
# Requirements
# ------------
# To run this tutorial, you need to install `aiida-workgraph` and `ase`. Open a terminal and run:
#
# .. code-block:: console
#
#    pip install aiida-workgraph ase
#
#
# Load the AiiDA profile.
#


from aiida import load_profile

load_profile()
#
# %%
# Atomization energy
# ===================================================
# Define a workgraph
# -------------------
# we use ASE emt calculator, but in practice you can use any calculator that is compatible with ASE.


# %%
# Second workflow: atomization energy molecule
# ==================================================
#
# The atomization energy, $\Delta E$, of a molecule can be expressed as:
#
# .. math::
#
#    \Delta E = n_{\text{atom}} \times E_{\text{atom}} - E_{\text{molecule}}
#
# Where:
#
# - :math:`\Delta E` is the atomization energy of the molecule.
# - :math:`n_{\text{atom}}` is the number of atoms.
# - :math:`E_{\text{atom}}` is the energy of an isolated atom.
# - :math:`E_{\text{molecule}}` is the energy of the molecule.
#
# Define a calcfunction to calculate the atomization energy
# ---------------------------------------------------------
#

from aiida_workgraph import task
from ase import Atoms
from ase.build import molecule


@task()
def calc_energy(atoms: Atoms) -> float:
    from ase.calculators.emt import EMT

    atoms.calc = EMT()
    atoms.get_potential_energy()
    return atoms.calc.results["energy"]


#
@task.graph()
def atomization_energy(molecule: Atoms, atoms: Atoms) -> float:
    """
    Calculate the atomization energy of a molecule.
    """
    e_atom = calc_energy(atoms).result
    e_molecule = calc_energy(molecule).result
    return len(molecule) * e_atom - e_molecule


# %%
# Submit the workgraph and print the atomization energy.
#
atoms = Atoms("N")
mol = molecule("N2")

wg = atomization_energy.build_graph(molecule=mol, atoms=atoms)
wg.to_html()

# %%
wg.run()
print(f"Atomization energy: {wg.outputs.result.value} eV")

# %%
# Let's have a look at the provenance graph.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Equations of state
# ==========================
# In this section, we will calculate the equation of state for a material using the `ase` package.
#


@task.graph(outputs=["energies", "volumes"])
def calc_all_energies(atoms: Atoms, strains: list) -> dict:
    """
    Calculate the equation of state for a material.
    """
    energies = {}
    volumes = {}

    for i, strain in enumerate(strains):
        strained_atoms = atoms.copy()
        strained_atoms.set_cell(strained_atoms.get_cell() * strain, scale_atoms=True)
        energies[f"strain_{i}"] = calc_energy(strained_atoms).result
        volumes[f"strain_{i}"] = strained_atoms.get_volume()

    return energies, volumes


@task()
def fit_eos(energies: dict, volumes: dict) -> dict:
    """Fit the EOS of the data."""
    from ase.eos import EquationOfState
    from ase.units import kJ

    volumes_list = list(volumes.values())
    energies_list = list(energies.values())
    eos = EquationOfState(volumes_list, energies_list)
    v0, e0, B = eos.fit()
    # convert B to GPa
    B = B / kJ * 1.0e24
    eos = {"energy unit": "eV", "v0": v0, "e0": e0, "B": B}
    return eos


@task.graph()
def eos_workflow(atoms: Atoms, strains: list) -> dict:
    """
    Calculate the equation of state for a material.
    """
    outputs = calc_all_energies(atoms, strains)
    eos = fit_eos(outputs.energies, outputs.volumes).result
    return eos


# %%
# Define the atoms and strains
from ase.build import bulk

cu = bulk("Cu", "fcc", a=3.6)
strains = [0.95, 0.98, 1.0, 1.02, 1.05]

# %%
# Build the workgraph and run it
wg_eos = eos_workflow.build_graph(atoms=cu, strains=strains)
wg_eos.to_html()

# %%
wg_eos.run()
print("Equation of state results: ", wg_eos.outputs.result.value)
