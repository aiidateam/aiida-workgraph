import pytest
from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain
from aiida_workgraph import task


def test_workchain(wg_workchain):
    """Submit simple calcjob."""
    wg = wg_workchain
    wg.name = "test_workchain"
    wg.run()
    # print("results: ", results[])
    assert wg.tasks.multiply_add2.outputs.result.value == 17


def test_build_workchain_inputs_outputs():
    """Submit simple calcjob."""

    node = task(MultiplyAddWorkChain)()._node
    inputs = MultiplyAddWorkChain.spec().inputs
    # inputs + metadata + _wait
    ninput = len(inputs.ports) + 1
    assert len(node.inputs) == ninput
    assert len(node.outputs) == 3


@pytest.mark.usefixtures("started_daemon_client")
def test_build_workchain(add_code):
    """Submit simple calcjob."""
    from aiida.orm import Int
    from aiida_workgraph import WorkGraph

    wg = WorkGraph(name="test_debug_math")
    code1 = wg.add_task("workgraph.load_code", "code1", pk=add_code.pk)
    multiply_add1 = wg.add_task(
        MultiplyAddWorkChain,
        "multiply_add1",
        x=Int(4),
        y=Int(2),
        z=Int(3),
    )
    multiply_add2 = wg.add_task(
        MultiplyAddWorkChain,
        "multiply_add2",
        x=Int(2),
        y=Int(3),
    )
    wg.add_link(code1.outputs[0], multiply_add1.inputs["code"])
    wg.add_link(code1.outputs[0], multiply_add2.inputs["code"])
    wg.add_link(multiply_add1.outputs[0], multiply_add2.inputs["z"])
    wg.submit(wait=True, timeout=100)
    assert wg.tasks.multiply_add2.outputs.result.value == 17
    # reload wg
    wg1 = WorkGraph.load(wg.pk)
    assert wg1.tasks.multiply_add2.outputs.result.value == 17
