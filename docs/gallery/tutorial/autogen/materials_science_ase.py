"""
=====================================================
Computational materials science with ASE
=====================================================

Introduction
============
We'll explore two key examples that highlight the flexibility of AiiDA-WorkGraph:

1.  **Atomization energy**: A simple, linear workflow to calculate the atomization energy of a diatomic molecule.
2.  **Equation of state (EOS)**: A more advanced workflow for bulk structure that showcases how to handle ``if`` condition, parallel execution, and dynamic inputs/outputs, a common pattern in dynamic workflows.

"""

# %%
# First, ensure you have `aiida-workgraph` and `ase` installed, and have configured an AiiDA environment.
# If not, you can install the necessary packages:
#
# .. code-block:: console
#
#    pip install aiida-workgraph ase
#
# Then, load your AiiDA profile.

from aiida import load_profile

load_profile()

# %%
# .. _ase_atomization_energy:
# Atomization energy of a diatomic molecule
# ==========================================
# The atomization energy (:math:`\Delta E`) is the energy required to break a molecule down into its individual, separate atoms.
# For a diatomic molecule, the formula is straightforward:
#
# .. math::
#
#    \Delta E = 2 \times E_{\text{atom}} - E_{\text{molecule}}
#
# Where :math:`E_{\text{atom}}` is the total energy of a single atom and :math:`E_{\text{molecule}}` is the total energy of the molecule.
#
# To build our workflow, we'll start with standard Python functions and convert them into AiiDA-trackable components using the ``@task`` decorator.
#
# .. note::
#
#    We're using the ASE EMT (Effective Medium Theory) calculator because it's exceptionally fast and perfect for demonstrations.
#    You can easily swap it with any other ASE-compatible calculator, like Quantum ESPRESSO, VASP, or GPAW, for your research.
#    For realistic simulations, especially with DFT codes, you would typically run these calculations on a remote computer.
#    For more details on running calculations remotely, please refer to the section on :ref:`Run calculations remotely <remote_calculations>`.
#
from aiida_workgraph import task, spec
from ase import Atoms
from ase.build import molecule


@task
def calculate_energy(atoms: Atoms) -> float:
    """Calculate the total energy of an atomic structure using ASE."""
    from ase.calculators.emt import EMT

    atoms.calc = EMT()
    atoms.get_potential_energy()
    return atoms.calc.results["energy"]


@task
def compute_atomization_energy(energy_atom: float, energy_molecule: float) -> float:
    """Calculate the atomization energy from atomic and molecular energies."""
    return 2 * energy_atom - energy_molecule


@task.graph()
def atomization_energy_workflow(molecule_obj: Atoms, atom_obj: Atoms) -> float:
    """Define the workflow graph to compute atomization energy."""

    e_atom = calculate_energy(atoms=atom_obj).result
    e_molecule = calculate_energy(atoms=molecule_obj).result
    return compute_atomization_energy(e_atom, e_molecule).result


# %%
# Build and run the workflow
# --------------------------


# %%
# First, we create the input structures for a nitrogen atom and molecule using ASE.
atom = Atoms("N")
mol = molecule("N2")

# %%
# Next, build the workgraph, but doesn't run it.
wg = atomization_energy_workflow.build(molecule_obj=mol, atom_obj=atom)

# %%
# You can visualize the planned workflow.
#
# .. note::
#
#    If you run in a Jupyter notebook, replace ``wg.to_html()`` with ``wg``.

wg.to_html()


# %%
# Now, execute the workgraph, which runs the tasks in the correct sequence.
wg.run()
print(f"Atomization energy for N2: {wg.outputs.result.value.value:.4f} eV")

# %%
# Visualize the Provenance Graph
# ------------------------------
# We can visualize the *provenance* graph of a completed workflow. This graph is the key to reproducibility,
# showing not just the tasks but also the actual data nodes that were created and stored in the AiiDA database.


wg.generate_provenance_graph()


# %%
# .. _ase_eos:
# Equation of state
# ==================
# Now for a more complex and practical example: calculating the Equation of State (EOS) for a bulk material.
# The process involves several steps:
#
# 1.  **Relax**(optional) the initial atomic structure to its lowest-energy state.
# 2.  **Strain** the relaxed structure by applying a series of scaling factors.
# 3.  **Calculate** the total energy and volume for each strained structure.
# 4.  **Fit** the resulting energy-volume data to an EOS model to find properties like the equilibrium volume and bulk modulus.
#
# This workflow perfectly demonstrates how to handle loops and a dynamic number of calculations.

from ase.calculators.emt import EMT
from ase.optimize import BFGS
from typing import Annotated


@task
def relax_structure(atoms: Atoms) -> Atoms:
    """Relax the atomic structure to its minimum energy configuration using ASE."""
    atoms.calc = EMT()
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.01)
    return atoms


# %%
# .. note::
#
#   If you want to run the ``relax_structure`` task on a remote computer, you can use the ``@task.pythonjob`` decorator.
#   Please refer to the section on :ref:`Run calculations remotely <remote_calculations>`.
#


@task
def create_strained_structures(
    atoms: Atoms, scales: list
) -> Annotated[dict, spec.namespace(scaled_structures=spec.dynamic(Atoms))]:
    """Generate a series of strained structures from a list of scaling factors."""
    scaled_structures = {}
    for i, scale in enumerate(scales):
        strained_atoms = atoms.copy()
        strained_atoms.set_cell(atoms.get_cell() * scale, scale_atoms=True)
        # Each structure gets a unique key, like "strain_0", "strain_1", etc.
        scaled_structures[f"strain_{i}"] = strained_atoms
    return {"scaled_structures": scaled_structures}


# %%
# .. note::
#
#    We emit each key in `scaled_structures` as a separate AiiDA output port by using a **dynamic namespace**.
#    For a full explanation of dynamic outputs and how to use namespaces, please refer to the section on :ref:`Dynamic Namespaces <dynamic_namespaces>`.


@task
def calculate_energy_and_volume(atoms: Atoms) -> dict:
    """Calculate the energy and volume for a single atomic structure."""
    atoms.calc = EMT()
    atoms.get_potential_energy()
    return {
        "energy": atoms.calc.results["energy"],
        "volume": atoms.get_volume(),
    }


@task.graph
def calc_all_structures(
    scaled_structures: Annotated[dict, spec.dynamic(Atoms)]
) -> Annotated[dict, spec.namespace(results=spec.dynamic(dict))]:
    """Sub-workflow to calculate energy and volume for all strained structures in parallel."""
    results = {}
    for key, atoms in scaled_structures.items():
        # The key for each result (e.g., "strain_0") becomes an output link
        # under the "results" namespace.
        results[key] = calculate_energy_and_volume(atoms).result

    # The returned dictionary's key "results" must match the name in the `outputs` decorator.
    return {"results": results}


@task
def fit_eos_model(data: Annotated[dict, spec.dynamic(dict)]) -> dict:
    """Fit Energy-Volume data to a Birch-Murnaghan Equation of State."""
    from ase.eos import EquationOfState
    from ase.units import kJ

    # Unpack the energies and volumes from the input data dictionary
    volumes_list = [value["volume"] for value in data.values()]
    energies_list = [value["energy"] for value in data.values()]

    eos = EquationOfState(volumes_list, energies_list)
    v0, e0, B = eos.fit()

    # The bulk modulus B is converted from eV/Å³ to GPa.
    B_GPa = B / kJ * 1.0e24
    return {"v0_A^3": v0, "e0_eV": e0, "B_GPa": B_GPa}


@task.graph()
def eos_workflow(atoms: Atoms, scales: list, run_relax: bool = True) -> dict:
    """The complete EOS workflow graph."""
    if run_relax:
        atoms = relax_structure(atoms=atoms).result
    strained = create_strained_structures(atoms=atoms, scales=scales)
    emt_outputs = calc_all_structures(scaled_structures=strained.scaled_structures)
    return fit_eos_model(data=emt_outputs.results).result


# %%
# Build and Run the EOS Workflow
# ------------------------------
# We first define the input crystal structure (fcc Copper) and the list of strains.
from ase.build import bulk

cu = bulk("Cu", "fcc", a=3.6)
scales = [0.95, 0.98, 1.0, 1.02, 1.05]

# %%
# Next, we build the workgraph with the inputs and visualize it:
wg = eos_workflow.build(atoms=cu, scales=scales)
wg.to_html()

# %%
# Finally, we run the workflow:

wg.run()

# %%
# The result is an AiiDA Dict node. We access its content via the `.value` attribute.
eos_result = wg.outputs.result.value
print("Equation of state results for Cu: ", eos_result.get_dict())

# %%
# Visualize the EOS Provenance Graph
# ----------------------------------
# This provenance graph is more complex, clearly showing the "fan-out" from
# ``create_strained_structures`` and the "fan-in" to ``fit_eos_model``, illustrating
# the power of AiiDA-WorkGraph to manage complex data flows automatically.


wg.generate_provenance_graph()

# %%
# Conclusion
# ==========
# Congratulations! You've now seen the core principles of `AiiDA-WorkGraph` in action.
# You have learned how to:
#
# -   Transform any Python function into a robust, provenance-tracked task with the ``@task`` decorator.
# -   Compose tasks into complete workflows using ``@task.graph``, from simple linear chains to complex graphs.
# -   Manage advanced patterns like ``if`` condition, parallel execution, and **dynamic namespaces**.
# -   Build, execute, and visualize both the workflow plan and its final, rich **provenance graph**.
#
# These powerful concepts are the foundation for building sophisticated, automated, and fully reproducible simulation pipelines for your own research projects. Happy computing! 🚀
