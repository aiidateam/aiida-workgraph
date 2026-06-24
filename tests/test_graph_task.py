from aiida_workgraph import task, WorkGraph, spec


@task
def add_multiply(a, b) -> spec.namespace(sum=int, product=int):
    """Task that returns a namespace with sum and product."""
    return {'sum': a + b, 'product': a * b}


def test_tuple_namespace_outputs():
    """Test mapping tuple results to the workgraph outputs."""

    @task.graph
    def test_graph(x, y) -> spec.namespace(
        out1=spec.namespace(sum=int, product=int),
        out2=spec.namespace(sum=int, product=int),
    ):
        return add_multiply(x, y), add_multiply(x, y)

    wg = test_graph.build(1, 2)
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

    results, wg = test_graph.run_get_graph(1, 2)
    #
    assert results['sum'] == 3
    assert results['product'] == 2
    # graph outputs
    assert wg.outputs.sum.value == 3
    assert wg.outputs.product.value == 2


@task
def make_label() -> str:
    return '777'


def test_scalar_input_behaves_as_plain_value():
    """A scalar input arrives in a ``@task.graph`` body as its plain Python value.

    Regression test for #786: at runtime the engine resolves the graph task's
    input to an ``orm.Str`` node. The body must be able to use it directly
    (``int(label)``) without reaching through ``.value``, mirroring how a plain
    Python scalar would behave.
    """

    @task.graph
    def inner(label) -> spec.namespace(sum=int, product=int):
        # `label` is delivered as orm.Str at runtime; int(label) must work
        # directly. Before the fix this raised TypeError.
        return add_multiply(int(label), 1)

    @task.graph
    def top() -> spec.namespace(sum=int, product=int):
        return inner(label=make_label().result)

    results, _ = top.run_get_graph()
    assert results['sum'] == 778
    assert results['product'] == 777


def test_scalar_input_forwarded_preserves_provenance():
    """Forwarding a scalar graph input unchanged links the original node.

    The body sees the raw value, but the input socket still holds the orm node,
    so forwarding the input to a sub-task must reuse the upstream node rather than
    serialise a fresh one. The unique value makes the node count unambiguous under
    the session-shared profile: exactly one such ``Str`` node proves the upstream
    node is reused, not duplicated. Regression guard for #786.
    """
    from aiida import orm

    unique = 'wg786-forward-provenance'

    @task
    def make_unique() -> str:
        return 'wg786-forward-provenance'

    @task
    def consume(x) -> int:
        return 1

    @task.graph
    def inner(label) -> int:
        # forward the scalar input unchanged (no .value); the link must carry the
        # upstream node, not a re-serialised copy of the raw value.
        return consume(x=label)

    @task.graph
    def top() -> int:
        return inner(label=make_unique().result)

    top.run_get_graph()
    query = orm.QueryBuilder().append(orm.Str, filters={'attributes.value': unique})
    assert query.count() == 1
