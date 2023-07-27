from aiida_worktree import build_node, WorkTree, node
from aiida import orm
from aiida.engine import calcfunction
from aiida.orm import Kind, Site, StructureData, KpointsData, load_node
from aiida import load_profile

load_profile()


@node()
@calcfunction
def get_marked_structures(structure, atoms_list, marker="X"):
    """"""
    site_info = [{"uuid": structure.uuid, "symbol": "ground"}]

    for data in atoms_list.get_list():
        index, orbital = data
        marked_structure = StructureData()
        kinds = {kind.name: kind for kind in structure.kinds}
        marked_structure.set_cell(structure.cell)

        for i, site in enumerate(structure.sites):
            if i == index:
                marked_kind = Kind(name=marker.value, symbols=site.kind_name)
                marked_site = Site(kind_name=marked_kind.name, position=site.position)
                marked_structure.append_kind(marked_kind)
                marked_structure.append_site(marked_site)
                symbol = site.kind_name
            else:
                if site.kind_name not in [kind.name for kind in marked_structure.kinds]:
                    marked_structure.append_kind(kinds[site.kind_name])
                new_site = Site(kind_name=site.kind_name, position=site.position)
                marked_structure.append_site(new_site)
        marked_structure.store()
        site_info.append(
            {
                "uuid": marked_structure.uuid,
                "index": index,
                "symbol": symbol,
                "orbital": orbital,
                "peak": f"{symbol}_{orbital}",
                "label": f"{symbol}_{index}",
            }
        )

    return orm.List(site_info)


# the structures is used to generate the worktree dynamically.
@node.group(outputs=[["gather1", "result", "result"]])
def run_scf(site_info, code, parameters, kpoints, pseudos, xps_pseudos, metadata):
    # register node
    ndata = {"path": "aiida_quantumespresso.calculations.pw.PwCalculation"}
    pw_node = build_node(ndata)
    #
    nt = WorkTree("run_scf")
    gather1 = nt.nodes.new("AiiDAGather", name="gather1")
    site_info = site_info.get_list()
    for site in site_info:
        element = site["symbol"]
        if element == "ground":
            continue
        pseudos[element] = xps_pseudos[f"{element}_gs"]
    # ground state
    structure_ground = load_node(site_info[0]["uuid"])
    pw_ground = nt.nodes.new(pw_node, name=f"ground")
    pw_ground.set(
        {
            "code": code,
            "parameters": parameters,
            "kpoints": kpoints,
            "pseudos": pseudos,
            "metadata": metadata,
            "structure": structure_ground,
        }
    )
    nt.links.new(pw_ground.outputs["output_parameters"], gather1.inputs[0])
    # excited state node
    for data in site_info[1:]:
        structure = load_node(data["uuid"])
        pseudos1 = pseudos.copy()
        pseudos1["X"] = xps_pseudos[data["peak"]]
        # remove pseudo of non-exist element
        pseudos1 = {kind.name: pseudos1[kind.name] for kind in structure.kinds}
        pw_excited = nt.nodes.new(pw_node, name=f"pw_excited_{data['label']}")
        pw_excited.set(
            {
                "code": code,
                "parameters": parameters,
                "kpoints": kpoints,
                "pseudos": pseudos1,
                "metadata": metadata,
                "structure": structure,
            }
        )
        nt.links.new(pw_excited.outputs["output_parameters"], gather1.inputs[0])
    return nt


# set link limit to a large value so that it can gather the result.
@node()
@calcfunction
def get_spectra(pw_outputs, site_info, correction_energies={}, orbital="1s"):
    from aiida.orm import load_node

    binding_energies = {}
    correction_energies = correction_energies.get_dict()
    ground = load_node(pw_outputs[0])
    n = len(pw_outputs)
    for i in range(1, n):
        # it only gather the uuid of the data, so we need to load it.
        data = load_node(pw_outputs[i])
        peak = site_info[i]["peak"]
        label = f'{site_info[i]["label"]}_{site_info[i]["orbital"]}'
        corr = correction_energies.get(peak, {})
        binding_energies[label] = (
            data.dict.energy
            - ground.dict.energy
            + corr.get("core", 0)
            - corr.get("exp", 0)
        )
    return orm.Dict(binding_energies)


@node.group(outputs=[["get_spectra1", "result", "result"]])
def xps(
    structure,
    code,
    atoms_list,
    parameters,
    kpoints,
    pseudos,
    xps_pseudos,
    metadata,
    correction_energies={},
):
    nt = WorkTree("xps")
    marked_structure1 = nt.nodes.new(
        get_marked_structures,
        name="marked_structure1",
        structure=structure,
        atoms_list=atoms_list,
    )
    run_scf1 = nt.nodes.new(run_scf, name="run_scf1")
    run_scf1.set(
        {
            "code": code,
            "parameters": parameters,
            "kpoints": kpoints,
            "pseudos": pseudos,
            "xps_pseudos": xps_pseudos,
            "metadata": metadata,
        }
    )
    get_spectra1 = nt.nodes.new(
        get_spectra, name="get_spectra1", correction_energies=correction_energies
    )
    nt.links.new(marked_structure1.outputs[0], run_scf1.inputs["site_info"])
    nt.links.new(run_scf1.outputs["result"], get_spectra1.inputs["pw_outputs"])
    nt.links.new(marked_structure1.outputs[0], get_spectra1.inputs["site_info"])
    return nt


# ===============================================================================
def load_xps_pseudo(pseudo_group="xps_pseudo_demo"):
    pseudo_group = (
        orm.QueryBuilder().append(orm.Group, filters={"label": pseudo_group}).one()[0]
    )
    pseudos = {node.label: node for node in pseudo_group.nodes}
    return pseudos, pseudo_group.base.extras.get("correction", {})


from ase.build import molecule
from aiida import orm

# create input structure node
mol = molecule("CH3CH2OH")
mol.pbc = True
mol.center(1.5)
mol = StructureData(ase=mol)
# create the PW node
code = orm.load_code("pw-7.2@localhost")
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
kpoints = KpointsData()
kpoints.set_kpoints_mesh([1, 1, 1])
# Load the pseudopotential family.
pseudo_family = orm.load_group("SSSP/1.2/PBE/efficiency")
pseudos = pseudo_family.get_pseudos(structure=mol)
#
xps_pseudos, correction_energies = load_xps_pseudo("xps_pseudo_demo")
#
metadata = {
    "options": {
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        },
    }
}
# ===============================================================================
nt = WorkTree("xps_test")
xps1 = nt.nodes.new(xps, name="xps")
xps1.set(
    {
        "structure": mol,
        "code": code,
        "atoms_list": [[0, "1s"], [2, "1s"]],
        "parameters": paras,
        "kpoints": kpoints,
        "pseudos": pseudos,
        "xps_pseudos": xps_pseudos,
        "metadata": metadata,
        "correction_energies": correction_energies,
    }
)
nt.submit(wait=True, timeout=300)
