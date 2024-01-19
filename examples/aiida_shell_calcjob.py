from aiida_worktree import node, WorkTree, build_node
from aiida import load_profile
from aiida.orm import load_code

load_profile()

pdb_fetch = load_code("pdb_fetch")
pdb_selchain = load_code("pdb_selchain")
pdb_delhetatm = load_code("pdb_delhetatm")
pdb_tidy = load_code("pdb_tidy")

shelljob = build_node({"path": "aiida_shell.calculations.shell.ShellJob"})


@node()
def generate_nodes(data):
    """Prepare the nodes"""
    return {"pdb": data}


# Create a worktree
wt = WorkTree(name="test_aiida_shell")
job1 = wt.nodes.new(shelljob, code=pdb_fetch, arguments=["1brs"])
job1.outputs.new("General", "stdout")
job2 = wt.nodes.new(shelljob, code=pdb_selchain, arguments=["-A,D", "{pdb}"])
job2.outputs.new("General", "stdout")
job3 = wt.nodes.new(shelljob, code=pdb_delhetatm, arguments=["{pdb}"])
job3.outputs.new("General", "stdout")
job4 = wt.nodes.new(shelljob, code=pdb_tidy, arguments=["{pdb}"])
job4.outputs.new("General", "stdout")
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
