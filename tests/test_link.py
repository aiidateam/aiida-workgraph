import aiida

aiida.load_profile()


def test_multiply_link():
    """Test multiply link."""

    from aiida_worktree import node, WorkTree
    from aiida.orm import Float, load_node
    from aiida.engine import calcfunction

    @node()
    @calcfunction
    def sum(inputs):
        total = 0
        for input in inputs:
            total += load_node(input).value
        return Float(total)

    nt = WorkTree(name="test_multiply_link")
    float1 = nt.nodes.new("AiiDANode", value=Float(1.0).store())
    float2 = nt.nodes.new("AiiDANode", value=Float(2.0).store())
    float3 = nt.nodes.new("AiiDANode", value=Float(3.0).store())
    gather1 = nt.nodes.new("AiiDAGather", "gather1")
    sum1 = nt.nodes.new(sum, "sum1")
    nt.links.new(float1.outputs[0], gather1.inputs[0])
    nt.links.new(float2.outputs[0], gather1.inputs[0])
    nt.links.new(float3.outputs[0], gather1.inputs[0])
    nt.links.new(gather1.outputs[0], sum1.inputs[0])
    nt.submit(wait=True)
    assert sum1.node.outputs.result.value == 6
