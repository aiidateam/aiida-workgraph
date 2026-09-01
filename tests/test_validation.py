from aiida_workgraph import task, namespace, meta, WorkGraph
from aiida_workgraph.errors import MissingRequiredInputsError
from aiida import orm
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


class TestMissingRequiredInputsPayload:
    """The missing-inputs failure carries structured entries next to the message."""

    def test_missing_leaf_input(self):
        @task.graph()
        def my_graph(a):
            add(a)

        with pytest.raises(MissingRequiredInputsError) as excinfo:
            my_graph.run(a=1)
        error = excinfo.value
        assert str(error).startswith('Missing required inputs:')
        assert 'add.y' in str(error)
        entries = {entry.socket_path: entry for entry in error.missing}
        assert set(entries) == {'add.y'}
        assert entries['add.y'].identifier == 'workgraph.any'
        assert entries['add.y'].help is None

    def test_missing_namespace_member(self):
        """A dotted namespace-member path resolves to its socket for identifier and help."""

        @task.graph()
        def my_graph(a, b: Annotated[dict, namespace(x=int, y=Annotated[int, meta(help='the second member')])]):
            add(a, b['x'])

        with pytest.raises(MissingRequiredInputsError) as excinfo:
            my_graph.run(a=1, b={'x': 1})
        entries = {entry.socket_path: entry for entry in excinfo.value.missing}
        assert set(entries) == {'graph_inputs.b.y'}
        assert entries['graph_inputs.b.y'].identifier == 'workgraph.int'
        assert entries['graph_inputs.b.y'].help == 'the second member'

    def test_help_text_is_carried(self):
        @task
        def add_documented(x, y: Annotated[int, meta(help='the right-hand addend')]):
            return x + y

        @task.graph()
        def my_graph(a):
            add_documented(x=a)

        with pytest.raises(MissingRequiredInputsError) as excinfo:
            my_graph.run(a=1)
        entries = {entry.socket_path: entry for entry in excinfo.value.missing}
        assert entries['add_documented.y'].help == 'the right-hand addend'

    def test_multiple_missing_inputs(self):
        @task.graph()
        def my_graph(a):
            add(a)
            add()

        with pytest.raises(MissingRequiredInputsError) as excinfo:
            my_graph.run(a=1)
        paths = [entry.socket_path for entry in excinfo.value.missing]
        assert paths == sorted(paths)
        assert set(paths) == {'add.y', 'add1.x', 'add1.y'}

    def test_error_is_a_value_error(self):
        assert issubclass(MissingRequiredInputsError, ValueError)

    def test_error_survives_pickle_and_deepcopy(self):
        import copy
        import pickle

        from aiida_workgraph.errors import MissingInput

        error = MissingRequiredInputsError([MissingInput('add.y', 'workgraph.int', 'The right operand.')])
        for clone in (pickle.loads(pickle.dumps(error)), copy.deepcopy(error)):
            assert clone.missing == error.missing
            assert str(clone) == str(error)

    def test_missing_code_input_identifier(self):
        @task
        def run_code(code: orm.InstalledCode, x):
            return x

        @task.graph()
        def my_graph(a):
            run_code(x=a)

        with pytest.raises(MissingRequiredInputsError) as excinfo:
            my_graph.run(a=1)
        entries = {entry.socket_path: entry for entry in excinfo.value.missing}
        assert entries['run_code.code'].identifier == 'workgraph.code'
