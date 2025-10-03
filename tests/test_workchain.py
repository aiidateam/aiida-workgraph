from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain
from node_graph.node_spec import SchemaSource
from aiida_workgraph import task


def test_build_workchain_inputs_outputs():
    """Run simple workchain."""

    node = task(MultiplyAddWorkChain)()._node
    inputs = MultiplyAddWorkChain.spec().inputs
    # inputs + metadata + _wait
    ninput = len(inputs.ports) + 1
    assert len(node.inputs) == ninput
    assert len(node.outputs) == 3
    assert node.spec.schema_source == SchemaSource.CALLABLE


def test_build_workchain(add_code):
    """Submit simple calcjob."""
    from aiida.orm import Int
    from aiida_workgraph import WorkGraph

    wg = WorkGraph(name='test_debug_math')
    wg.add_task(
        MultiplyAddWorkChain,
        'multiply_add1',
        x=Int(4),
        y=Int(2),
        z=Int(3),
        code=add_code,
    )
    wg.run()
    assert wg.tasks.multiply_add1.outputs.result.value == 11
    # reload wg
    wg1 = WorkGraph.load(wg.pk)
    assert wg1.tasks.multiply_add1.get_executor().callable == MultiplyAddWorkChain
