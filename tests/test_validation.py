from aiida_workgraph import task, namespace, WorkGraph
from typing import Annotated
import pytest
import re


@task
def add(x, y):
    return x + y


def test_validate_required_inputs():
    @task.graph()
    def my_graph(a, b: Annotated[dict, namespace(x=int, y=int)]):
        add(a, b['x'])
        add(a)

    with pytest.raises(
        ValueError,
        match=re.escape('Missing required inputs:'),
    ):
        my_graph.run(a=1, b={'x': 1})


# The default makes the namespace optional; it is deliberately not ``None``.
# Up to Python 3.10 ``get_type_hints`` rewrites a ``None``-defaulted annotation
# as ``Optional[...]``, and deduplicating that union hashes the ``Annotated``
# metadata — which holds a dict, so the spec fails to build.
@task
def consume_optional_namespace(payload: Annotated[dict, namespace(x=int, y=int)] = {}):
    return payload


@task
def consume_required_namespace(payload: Annotated[dict, namespace(x=int, y=int)]):
    return payload


def test_optional_namespace_left_empty_is_not_missing():
    """An optional namespace socket that the caller never touched is complete.

    Its children are required only *within* the namespace — they say what must be
    supplied if the namespace is supplied at all. Recursing into an untouched
    optional namespace reported every child as missing, which made an optional
    grouped input (e.g. a ``TypedDict | None``) impossible to omit.
    """
    wg = WorkGraph('optional-namespace')
    wg.add_task(consume_optional_namespace, name='consumer')

    assert wg.find_missing_inputs(wg.tasks.consumer.inputs) == []
    wg.check_required_inputs()


def test_partially_filled_optional_namespace_reports_its_gaps():
    """Once an optional namespace is used at all, its children are back in force."""
    wg = WorkGraph('optional-namespace-partial')
    task_ = wg.add_task(consume_optional_namespace, name='consumer')
    task_.inputs.payload.x.value = 1

    assert wg.find_missing_inputs(task_.inputs) == ['consumer.payload.y']
    with pytest.raises(ValueError, match=re.escape('consumer.payload.y')):
        wg.check_required_inputs()


@task
def produce_namespace() -> Annotated[dict, namespace(x=int, y=int)]:
    return {'x': 1, 'y': 2}


def test_optional_namespace_with_a_linked_child_reports_its_gaps():
    """A link into a single child counts as using the namespace.

    The link sits on the child, so the namespace itself has neither a value nor a
    link of its own; deciding on those alone would skip it and hide the sibling
    that the caller still has to provide.
    """
    wg = WorkGraph('optional-namespace-linked-child')
    producer = wg.add_task(produce_namespace, name='producer')
    consumer = wg.add_task(consume_optional_namespace, name='consumer')
    consumer.inputs.payload.x = producer.outputs.x

    assert wg.find_missing_inputs(consumer.inputs) == ['consumer.payload.y']


def test_required_namespace_left_empty_still_reports_its_children():
    """A required namespace is unaffected: omitting it is still an error."""
    wg = WorkGraph('required-namespace')
    task_ = wg.add_task(consume_required_namespace, name='consumer')

    assert sorted(wg.find_missing_inputs(task_.inputs)) == [
        'consumer.payload.x',
        'consumer.payload.y',
    ]


@pytest.mark.parametrize(
    'name, reason',
    [
        pytest.param('_hidden', 'cannot start with an underscore', id='leading-underscore'),
        pytest.param('hidden_', 'cannot end with an underscore', id='trailing-underscore'),
        pytest.param('2nd', 'not a valid python identifier', id='non-identifier'),
    ],
)
def test_invalid_task_name_raises_at_build_time(name, reason):
    """An explicit, invalid ``name=`` to the low-level ``add_task`` is rejected at build time.

    Regression test for https://github.com/aiidateam/aiida-workgraph/issues/784: such names
    previously passed the (weaker) node_graph name check and only failed at run time inside
    the engine, where the failure was swallowed.
    """
    wg = WorkGraph()
    with pytest.raises(ValueError, match=re.escape(f"Invalid task name '{name}'")) as excinfo:
        wg.add_task(add, name=name)
    message = str(excinfo.value)
    assert reason in message
    # the fix hint is the one for the low-level API, not the call_link_label one
    assert 'WorkGraph.add_task' in message


def test_invalid_derived_task_name_raises_at_build_time():
    """The issue #784 example: an invalid name derived from the function name (never passed
    explicitly), built through the high-level ``@task.graph`` API, must also be caught at
    build time. This path still resolves the name inside ``add_task``.
    """

    @task
    def _hidden() -> dict:
        return {'ran': True}

    @task.graph
    def top():
        _hidden()

    with pytest.raises(ValueError, match=re.escape("Invalid task name '_hidden'")) as excinfo:
        top.build()
    message = str(excinfo.value)
    # the name was derived, so the fix is to rename the callable; the user never called
    # `add_task` here, so the message must not tell them to use it
    assert 'rename the function/callable' in message
    assert 'WorkGraph.add_task' not in message


def test_invalid_call_link_label_raises_at_build_time():
    """An invalid name set via ``metadata={'call_link_label': ...}`` must be rejected too.

    This is the documented way to give a task an explicit name in the high-level API. The
    override is applied in ``TaskHandle.__call__`` *after* ``add_task`` validated the derived
    name, so it needs its own check; without it the invalid label slips through and only
    fails silently at run time (issue #784). The user-facing fix hint here is the
    ``call_link_label`` one, not the low-level ``name=`` one.
    """
    with pytest.raises(ValueError, match=re.escape("Invalid task name '_sum'")) as excinfo:
        with WorkGraph():
            add(1, 2, metadata={'call_link_label': '_sum'})
    assert 'call_link_label' in str(excinfo.value)
