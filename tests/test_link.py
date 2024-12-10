import pytest


@pytest.mark.usefixtures("started_daemon_client")
def test_multiply_link() -> None:
    """Test multiply link."""

    from aiida_workgraph import task, WorkGraph
    from aiida.orm import Float

    @task.calcfunction()
    def sum(**datas):
        total = 0
        for _, data in datas.items():
            total += data.value
        return Float(total)

    wg = WorkGraph(name="test_multiply_link")
    float1 = wg.add_task("workgraph.aiida_node", pk=Float(1.0).store().pk)
    float2 = wg.add_task("workgraph.aiida_node", pk=Float(2.0).store().pk)
    float3 = wg.add_task("workgraph.aiida_node", pk=Float(3.0).store().pk)
    sum1 = wg.add_task(sum, "sum1")
    sum1.inputs["datas"]._link_limit = 100
    wg.add_link(float1.outputs[0], sum1.inputs["datas"])
    wg.add_link(float2.outputs[0], sum1.inputs["datas"])
    wg.add_link(float3.outputs[0], sum1.inputs["datas"])
    # wg.submit(wait=True)
    wg.run()
    assert sum1.node.outputs.result.value == 6
