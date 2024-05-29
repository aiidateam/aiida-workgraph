from aiida_workgraph import node, WorkGraph
from aiida import load_profile
from aiida.orm import load_code

load_profile()

pdb_fetch = load_code("pdb_fetch")
pdb_selchain = load_code("pdb_selchain")
pdb_delhetatm = load_code("pdb_delhetatm")
pdb_tidy = load_code("pdb_tidy")


@node()
def generate_nodes(file):
    """Prepare the nodes"""
    return {"pdb": file}


# Create a workgraph
wg = WorkGraph(name="test_aiida_shell_calcjob")
job1 = wg.nodes.new("ShellJob", code=pdb_fetch, arguments=["1brs"])
job2 = wg.nodes.new("ShellJob", code=pdb_selchain, arguments=["-A,D", "{pdb}"])
job3 = wg.nodes.new("ShellJob", code=pdb_delhetatm, arguments=["{pdb}"])
job4 = wg.nodes.new("ShellJob", code=pdb_tidy, arguments=["{pdb}"])
generate_nodes1 = wg.nodes.new(generate_nodes)
generate_nodes2 = wg.nodes.new(generate_nodes)
generate_nodes3 = wg.nodes.new(generate_nodes)
wg.links.new(job1.outputs["stdout"], generate_nodes1.inputs[0])
wg.links.new(generate_nodes1.outputs[0], job2.inputs["nodes"])
wg.links.new(job2.outputs["stdout"], generate_nodes2.inputs[0])
wg.links.new(generate_nodes2.outputs[0], job3.inputs["nodes"])
wg.links.new(job3.outputs["stdout"], generate_nodes3.inputs[0])
wg.links.new(generate_nodes3.outputs[0], job4.inputs["nodes"])
wg.submit()
