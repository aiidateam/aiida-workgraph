from aiida_workgraph import task, WorkGraph


def test_tuple_outputs(decorated_add, decorated_multiply):
    """Test mapping tuple results to the workgraph outputs."""

    @task.graph(outputs=["sum", "product"])
    def add_multiply(a, b):
        """Task that returns the sum and product of two numbers."""
        sum_result = decorated_add(a, b).result
        product_result = decorated_multiply(sum_result, b).result
        return sum_result, product_result

    wg = add_multiply.build_graph(1, 2)
    assert "sum" in wg.outputs
    assert "product" in wg.outputs

    with WorkGraph() as wg:
        outputs = add_multiply(1, 2)
        wg.run()
    assert outputs.sum.value == 3
    assert outputs.product.value == 6
