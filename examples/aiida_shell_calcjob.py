from aiida_worktree import node, WorkTree
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


# Create a worktree
wt = WorkTree(name="test_aiida_shell_calcjob")
job1 = wt.nodes.new("AiiDAShell", code=pdb_fetch, arguments=["1brs"])
job2 = wt.nodes.new("AiiDAShell", code=pdb_selchain, arguments=["-A,D", "{pdb}"])
job3 = wt.nodes.new("AiiDAShell", code=pdb_delhetatm, arguments=["{pdb}"])
job4 = wt.nodes.new("AiiDAShell", code=pdb_tidy, arguments=["{pdb}"])
generate_nodes1 = wt.nodes.new(generate_nodes)
generate_nodes2 = wt.nodes.new(generate_nodes)
generate_nodes3 = wt.nodes.new(generate_nodes)
wt.links.new(job1.outputs["stdout"], generate_nodes1.inputs[0])
wt.links.new(generate_nodes1.outputs[0], job2.inputs["nodes"])
wt.links.new(job2.outputs["stdout"], generate_nodes2.inputs[0])
wt.links.new(generate_nodes2.outputs[0], job3.inputs["nodes"])
wt.links.new(job3.outputs["stdout"], generate_nodes3.inputs[0])
wt.links.new(generate_nodes3.outputs[0], job4.inputs["nodes"])
wt.submit()
