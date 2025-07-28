"""
=====================================================
Computational materials science with ASE
=====================================================

Introduction
============
Unlock the power of automated and reproducible computational science!
This tutorial will guide you through building, running, and visualizing computational workflows using `AiiDA-WorkGraph` and the `Atomistic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_.
You'll learn how to construct complex pipelines where every calculation and data point is automatically tracked, ensuring your research is completely reproducible.

We'll explore two key examples that highlight the flexibility of AiiDA-WorkGraph:

1.  **Atomization energy**: A simple, linear workflow to calculate the atomization energy of a diatomic molecule.
2.  **Equation of state (EOS)**: A more advanced workflow for bulk structure that showcases how to handle loops and dynamic inputs/outputs, a common pattern in dynamic workflows.

"""

# %%
# First, ensure you have `aiida-workgraph`, `ase`, and a configured AiiDA environment.
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
# Atomization energy of a diatomic molecule
# ==========================================
# The atomization energy (ΔE) is the energy required to break a molecule down into its individual, separate atoms.
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
#
from aiida_workgraph import task
from ase import Atoms
from ase.build import molecule


@task()
def calculate_energy(atoms: Atoms) -> float:
    """Calculate the total energy of an atomic structure using ASE."""
    from ase.calculators.emt import EMT

    atoms.calc = EMT()
    atoms.get_potential_energy()
    return atoms.calc.results["energy"]


@task()
def compute_atomization_energy(energy_atom: float, energy_molecule: float) -> float:
    """Calculate the atomization energy from atomic and molecular energies."""
    result = 2 * energy_atom - energy_molecule
    return result


# %%
# .. note::
#
#    When you call a task function like ``calculate_energy(...)`` inside a workgraph, it **doesn't run immediately**.
#    Instead, it creates a task node in the workflow and returns a *reference* of its future result.
#    You can then wire this reference as an input to the next task, defining the dependencies between tasks.


@task.graph()
def atomization_energy_workflow(molecule_obj: Atoms, atom_obj: Atoms) -> float:
    """Define the workflow graph to compute atomization energy."""

    e_atom = calculate_energy(atoms=atom_obj).result
    e_molecule = calculate_energy(atoms=molecule_obj).result
    result = compute_atomization_energy(e_atom, e_molecule).result
    return result


# %%
# Build and Run the Workflow
# --------------------------
# First, we create the input structures for a nitrogen atom and molecule using ASE.
atom = Atoms("N")
mol = molecule("N2")

# Next, build the workgraph, but doesn't run it.
wg = atomization_energy_workflow.build_graph(molecule_obj=mol, atom_obj=atom)

# You can generate an HTML file to visualize the planned workflow.
wg.to_html()
# If you're in a Jupyter notebook, you can often display the graph directly.
# wg

# %%
# Now, execute the workgraph, which runs the tasks in the correct sequence.
wg.run()
print(f"Atomization energy for N2: {wg.outputs.result.value.value:.4f} eV")

# %%
# Visualize the Provenance Graph
# ------------------------------
# After the run, we can visualize the *provenance* graph. This graph is the key to reproducibility,
# showing not just the tasks but also the actual data nodes that were created and stored in the AiiDA database.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)


# %%
# Equation of state
# ==================
# Now for a more complex and practical example: calculating the Equation of State (EOS) for a bulk material.
# The process involves several steps:
#
# 1.  **Relax** the initial atomic structure to its lowest-energy state.
# 2.  **Strain** the relaxed structure by applying a series of scaling factors.
# 3.  **Calculate** the total energy and volume for each strained structure.
# 4.  **Fit** the resulting energy-volume data to an EOS model to find properties like the equilibrium volume and bulk modulus.
#
# This workflow perfectly demonstrates how to handle loops and a dynamic number of calculations.

from ase.calculators.emt import EMT
from ase.optimize import BFGS


@task()
def relax_structure(atoms: Atoms) -> Atoms:
    """Relax the atomic structure to its minimum energy configuration using ASE."""
    atoms.calc = EMT()
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.01)
    return atoms


@task(
    outputs={
        "scaled_structures": {
            "identifier": "workgraph.namespace",
            "metadata": {"dynamic": True},
        }
    }
)
def create_strained_structures(atoms: Atoms, scales: list) -> dict:
    """Generate a series of strained structures from a list of scaling factors."""
    scaled_structures = {}
    for i, scale in enumerate(scales):
        strained_atoms = atoms.copy()
        strained_atoms.set_cell(atoms.get_cell() * scale, scale_atoms=True)
        # Each structure gets a unique key, like "strain_0", "strain_1", etc.
        scaled_structures[f"strain_{i}"] = strained_atoms
    return {"scaled_structures": scaled_structures}


# %%
# .. important::
#    **Handling dynamic outputs with namespaces** ✨
#
#    In our EOS workflow, the ``create_strained_structures`` task generates a *variable* number of outputs. How do we track each one?
#
#    The solution is a **dynamic output namespace**. By adding the ``outputs={...}`` argument to the ``@task`` decorator, we tell AiiDA-WorkGraph:
#    "The ``scaled_structures`` output is a dictionary, but I want you to treat *each key-value pair* as a separate, named output port (e.g., ``scaled_structures.strain_0``, ``scaled_structures.strain_1``)."
#
#    This is crucial for two reasons:
#
#    1.  **Full Provenance**: Every strained structure becomes an individual, tracked AiiDA node.
#        This gives you a crystal-clear history of your entire calculation, with no "black boxes".
#    2.  **Database Integrity**: It avoids saving a list of ``Atoms`` objects as a single, opaque file (e.g., using pickle, disabled by default).
#        Instead, each structure is a queryable, first-class citizen in the AiiDA database.


@task()
def calculate_energy_and_volume(atoms: Atoms) -> dict:
    """Calculate the energy and volume for a single atomic structure."""
    atoms.calc = EMT()
    atoms.get_potential_energy()
    results = {
        "energy": atoms.calc.results["energy"],
        "volume": atoms.get_volume(),
    }
    return results


@task.graph(
    outputs={
        "results": {"identifier": "workgraph.namespace", "metadata": {"dynamic": True}}
    }
)
def calc_all_structures(**scaled_structures) -> dict:
    """Sub-workflow to calculate energy and volume for all strained structures in parallel."""
    results = {}
    for key, atoms in scaled_structures.items():
        # The key for each result (e.g., "strain_0") becomes an output link
        # under the "results" namespace.
        results[key] = calculate_energy_and_volume(atoms).result

    # The returned dictionary's key "results" must match the name in the `outputs` decorator.
    return {"results": results}


# %%
# .. note::
#
#    The ``**scaled_structures`` syntax in the ``calc_all_structures`` function signature is Python's way of accepting an arbitrary number of keyword arguments.
#    This makes it the perfect receiver for the dynamic outputs from ``create_strained_structures``, elegantly handling the ``strain_0=...``, ``strain_1=...``, etc., inputs.


@task()
def fit_eos_model(**data) -> dict:
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
    result = {"v0_A^3": v0, "e0_eV": e0, "B_GPa": B_GPa}
    return result


@task.graph()
def eos_workflow(atoms: Atoms, scales: list) -> dict:
    """The complete EOS workflow graph."""
    relaxed_atoms = relax_structure(atoms=atoms).result
    strained = create_strained_structures(atoms=relaxed_atoms, scales=scales)
    emt_outputs = calc_all_structures(scaled_structures=strained.scaled_structures)
    result = fit_eos_model(data=emt_outputs.results).result
    return result


# %%
#
# .. important::
#
#    When linking task outputs to a keyword argument of the task, you must pass the **output sockets** explicitly by name, e.g.:
#
#    .. code-block:: python
#
#        calc_all_structures(scaled_structures=strained.scaled_structures)
#        fit_eos_model(data=emt_outputs.results).result
#
#    rather than using argument unpacking like:
#
#    .. code-block:: python
#
#        calc_all_structures(**strained.scaled_structures)
#        fit_eos_model(**emt_outputs.results)
#
#    Because ``strained.scaled_structures`` and ``emt_outputs.results`` are *references* to the future output dictionaries (i.e. output sockets), not the actual values.
#    Naming the arguments ensures AiiDA-WorkGraph knows exactly which output port to connect to which input port.
#

# %%
# Build and Run the EOS Workflow
# ------------------------------
# Define the input crystal structure (fcc Copper) and the list of strains.
from ase.build import bulk

cu = bulk("Cu", "fcc", a=3.6)
scales = [0.95, 0.98, 1.0, 1.02, 1.05]

# %%
# Build the workgraph and visualize the plan.
wg = eos_workflow.build_graph(atoms=cu, scales=scales)
wg.to_html()

# %%
# Run the workflow.
wg.run()

# The result is an AiiDA Dict node. We access its content via the `.value` attribute.
eos_result = wg.outputs.result.value
print("Equation of state results for Cu: ", eos_result.value)

# %%
# Visualize the EOS Provenance Graph
# ----------------------------------
# This provenance graph is more complex, clearly showing the "fan-out" from
# `create_strained_structures` and the "fan-in" to `fit_eos_model`, illustrating
# the power of AiiDA-WorkGraph to manage complex data flows automatically.

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Conclusion
# ==========
# Congratulations! You've now seen the core principles of `AiiDA-WorkGraph` in action.
# You have learned how to:
#
# -   Transform any Python function into a robust, provenance-tracked task with the ``@task`` decorator.
# -   Compose tasks into complete workflows using ``@task.graph``, from simple linear chains to complex graphs.
# -   Manage advanced patterns like loops and parallel execution using **dynamic namespaces**.
# -   Build, execute, and visualize both the workflow plan and its final, rich **provenance graph**.
#
# These powerful concepts are the foundation for building sophisticated, automated, and fully reproducible simulation pipelines for your own research projects. Happy computing! 🚀
