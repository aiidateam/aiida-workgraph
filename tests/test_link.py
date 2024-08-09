import pytest


@pytest.mark.usefixtures("started_daemon_client")
def test_multiply_link() -> None:
    """Test multiply link."""

    from aiida_workgraph import task, WorkGraph
    from aiida.orm import Float, load_node

    @task.calcfunction()
    def sum(datas):
        total = 0
        for data in datas:
            total += load_node(data).value
        return Float(total)

    wg = WorkGraph(name="test_multiply_link")
    float1 = wg.add_task("workgraph.aiida_node", pk=Float(1.0).store().pk)
    float2 = wg.add_task("workgraph.aiida_node", pk=Float(2.0).store().pk)
    float3 = wg.add_task("workgraph.aiida_node", pk=Float(3.0).store().pk)
    gather1 = wg.add_task("workgraph.gather", "gather1")
    sum1 = wg.add_task(sum, "sum1")
    wg.add_link(float1.outputs[0], gather1.inputs[0])
    wg.add_link(float2.outputs[0], gather1.inputs[0])
    wg.add_link(float3.outputs[0], gather1.inputs[0])
    wg.add_link(gather1.outputs["result"], sum1.inputs["datas"])
    wg.submit(wait=True)
    assert sum1.node.outputs.result.value == 6
