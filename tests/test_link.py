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

    wt = WorkTree(name="test_multiply_link")
    float1 = wt.nodes.new("AiiDANode", value=Float(1.0).store())
    float2 = wt.nodes.new("AiiDANode", value=Float(2.0).store())
    float3 = wt.nodes.new("AiiDANode", value=Float(3.0).store())
    gather1 = wt.nodes.new("AiiDAGather", "gather1")
    sum1 = wt.nodes.new(sum, "sum1")
    wt.links.new(float1.outputs[0], gather1.inputs[0])
    wt.links.new(float2.outputs[0], gather1.inputs[0])
    wt.links.new(float3.outputs[0], gather1.inputs[0])
    wt.links.new(gather1.outputs[0], sum1.inputs[0])
    wt.submit(wait=True)
    assert sum1.node.outputs.result.value == 6
