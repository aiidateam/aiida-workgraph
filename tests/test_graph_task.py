from aiida_workgraph import task, WorkGraph, spec


@task
def add_multiply(a, b) -> spec.namespace(sum=int, product=int):
    """Task that returns a namespace with sum and product."""
    return {"sum": a + b, "product": a * b}


def test_tuple_namespace_outputs():
    """Test mapping tuple results to the workgraph outputs."""

    @task.graph
    def test_graph(
        x, y
    ) -> spec.namespace(
        out1=spec.namespace(sum=int, product=int),
        out2=spec.namespace(sum=int, product=int),
    ):
        return add_multiply(x, y), add_multiply(x, y)

    wg = test_graph.build_graph(1, 2)
    # graph inputs
    assert wg.inputs.x.value == 1
    assert wg.inputs.y.value == 2
    wg.run()
    # graph outputs
    assert wg.outputs.out1.sum.value == 3
    assert wg.outputs.out1.product.value == 2

    with WorkGraph() as wg:
        outputs = test_graph(1, 2)
        wg.run()
    assert outputs.out1.sum.value == 3
    assert outputs.out2.product.value == 2


def test_single_namespace_outputs():
    """Test mapping namespace results to the workgraph outputs."""

    @task.graph
    def test_graph(x, y) -> spec.namespace(sum=int, product=int):
        return add_multiply(x, y)

    wg = test_graph.build_graph(1, 2)
    wg.run()
    # graph outputs
    assert wg.outputs.sum.value == 3
    assert wg.outputs.product.value == 2
