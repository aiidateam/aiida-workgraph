import pytest
from aiida_workgraph import WorkGraph, task, spec
from aiida_workgraph.workgraph import BOOKKEEPING_KEYS, ENGINE_LAUNCH_KEYS, WorkGraphMetadata
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from typing import Any
import re


def test_represent():
    """Test the __repr__ method of WorkGraph."""
    wg = WorkGraph('test_represent')
    assert repr(wg) == f'WorkGraph(name="test_represent", uuid="{wg.uuid}")'
    assert str(wg) == f'WorkGraph(name="test_represent", uuid="{wg.uuid}")'


def test_should_submit():
    wg = WorkGraph('test_should_submit')
    wg.save()
    with pytest.raises(
        ValueError, match=re.escape(f'Process {wg.pk} has already been created. Please use the submit() method.')
    ):
        wg.run()


def test_from_dict(decorated_add):
    """Export NodeGraph to dict."""
    wg = WorkGraph('test_from_dict')
    task1 = wg.add_task(decorated_add, x=2, y=3)
    wg.add_task('workgraph.test_sum_diff', name='sumdiff2', x=4, y=task1.outputs.result)
    wgdata = wg.to_dict()
    wg1 = WorkGraph.from_dict(wgdata)
    assert len(wg.tasks) == len(wg1.tasks)
    assert len(wg.links) == len(wg1.links)


def test_add_task():
    """Add add task."""
    wg = WorkGraph('test_add_task')
    add1 = wg.add_task(ArithmeticAddCalculation, name='add1')
    add2 = wg.add_task(ArithmeticAddCalculation, name='add2')
    wg.add_link(add1.outputs.sum, add2.inputs.x)
    assert len(wg.tasks) == 5
    assert len(wg.links) == 1


def test_show_state(wg_task):
    from io import StringIO
    import sys

    # Redirect stdout to capture prints
    captured_output = StringIO()
    sys.stdout = captured_output
    # Call the method
    wg_task.name = 'test_show_state'
    wg_task.show()
    # Reset stdout
    sys.stdout = sys.__stdout__
    # Check the output
    output = captured_output.getvalue()
    assert 'WorkGraph: test_show_state, PK: None, State: CREATED' in output
    assert 'sumdiff1' in output
    assert 'PLANNED' in output


def test_save_load(wg_task, decorated_add):
    """Save the workgraph"""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
    from aiida_workgraph.executors.builtins import UnavailableExecutor

    wg: WorkGraph = wg_task
    wg.add_task(decorated_add, name='add1', x=2, y=3)
    metadata = {
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 2,
            },
        }
    }
    wg.add_task(
        ArithmeticAddCalculation,
        name='add2',
        x=4,
        y=wg.tasks.add1.outputs.result,
        metadata=metadata,
    )
    wg.name = 'test_save_load'
    wg.save()
    assert wg.process.process_state.value.upper() == 'CREATED'
    assert wg.process.process_label == 'WorkGraph<test_save_load>'
    assert wg.process.label == 'test_save_load'

    wg2 = WorkGraph.load(wg.process.pk)
    assert len(wg.tasks) == len(wg2.tasks)
    # the executor of the decorated task is pickled,
    # so it's not stored in the database.
    assert wg2.tasks.add1.get_executor().callable == UnavailableExecutor
    # The ArithmeticAddCalculation is importable,
    # so we can restore the executor from the module path.
    assert wg2.tasks.add2.get_executor().callable == wg.tasks.add2.get_executor().callable
    assert wg.tasks.add2.inputs.metadata._value == wg2.tasks.add2.inputs.metadata._value
    # metadata is also loaded
    assert wg2.tasks.add2.inputs.metadata.options.resources.value['num_mpiprocs_per_machine'] == 2
    # TODO, the following code is not working
    # wg2.save()
    # assert wg2.tasks.add1.executor == decorated_add
    # remove the extra, will raise an error

    # Check that it also works for the uuid
    wg3 = WorkGraph.load(wg.process.uuid)
    assert len(wg.tasks) == len(wg3.tasks)

    assert wg3.tasks.add1.get_executor().callable == UnavailableExecutor

    assert wg3.tasks.add2.get_executor().callable == wg.tasks.add2.get_executor().callable
    assert wg.tasks.add2.inputs.metadata._value == wg3.tasks.add2.inputs.metadata._value

    assert wg3.tasks.add2.inputs.metadata.options.resources.value['num_mpiprocs_per_machine'] == 2


def test_explicit_label_survives_save(wg_task):
    """An explicit ``metadata.label`` must not be overwritten by the workgraph name."""
    wg = wg_task
    wg.name = 'test_explicit_label_survives_save'
    wg.save(metadata={'label': 'my-explicit-label'})
    assert wg.process.process_label == 'WorkGraph<test_explicit_label_survives_save>'
    assert wg.process.label == 'my-explicit-label'


def test_empty_label_survives_save(wg_task):
    """An explicit empty-string ``metadata.label`` is a given label, not an absent one.

    aiida-core's ``Process._setup_metadata`` applies ``label`` on key presence, not
    truthiness, so ``metadata={'label': ''}`` must not fall back to the workgraph name.
    """
    wg = wg_task
    wg.name = 'test_empty_label_survives_save'
    wg.save(metadata={'label': ''})
    assert wg.process.label == ''


def test_explicit_label_survives_run(wg_task):
    """An explicit ``metadata.label`` passed to ``run()`` must survive on the process node."""
    wg = wg_task
    wg.name = 'test_explicit_label_survives_run'
    wg.run(metadata={'label': 'my-explicit-label'})
    assert wg.process.process_label == 'WorkGraph<test_explicit_label_survives_run>'
    assert wg.process.label == 'my-explicit-label'


@pytest.mark.usefixtures('started_daemon_client')
def test_explicit_label_survives_submit(wg_task):
    """An explicit ``metadata.label`` passed to ``submit()`` must survive on the process node."""
    wg = wg_task
    wg.name = 'test_explicit_label_survives_submit'
    wg.submit(metadata={'label': 'my-explicit-label'}, wait=True, timeout=30)
    assert wg.process.process_label == 'WorkGraph<test_explicit_label_survives_submit>'
    assert wg.process.label == 'my-explicit-label'


def test_explicit_label_survives_nested_graph():
    """An explicit ``metadata.label`` on a nested ``@task.graph`` call survives to the child node."""

    @task()
    def add(x, y):
        return x + y

    @task.graph()
    def inner_graph(x, y):
        return add(x, y).result

    with WorkGraph('test_explicit_label_survives_nested_graph') as wg:
        inner_graph(1, 2, metadata={'label': 'nested-explicit-label'})
        wg.run()

    called = wg.process.called
    assert len(called) == 1
    # the child's process_label is derived from the call-site task name (here the
    # function's default name), not from the explicit label just asserted below.
    assert called[0].process_label == 'WorkGraph<inner_graph>'
    assert called[0].label == 'nested-explicit-label'


def test_wg_metadata_forwarded_to_run_via_constructor():
    """A ``label`` set through the constructor's ``metadata`` kwarg reaches ``process.label`` on ``run()``."""
    wg = WorkGraph('test_wg_metadata_forwarded_to_run_via_constructor', metadata={'label': 'Human-readable label'})
    wg.run()
    assert wg.process.process_label == 'WorkGraph<test_wg_metadata_forwarded_to_run_via_constructor>'
    assert wg.process.label == 'Human-readable label'


def test_wg_metadata_forwarded_to_run_via_attribute(wg_task):
    """A ``label`` set through the ``wg.metadata`` attribute reaches ``process.label`` on ``run()``."""
    wg = wg_task
    wg.name = 'test_wg_metadata_forwarded_to_run_via_attribute'
    wg.metadata['label'] = 'Human-readable label'
    wg.run()
    assert wg.process.label == 'Human-readable label'


@pytest.mark.usefixtures('started_daemon_client')
def test_wg_metadata_forwarded_to_submit(wg_task):
    """``wg.metadata['label']`` is forwarded to the process node's ``label`` on ``submit()``."""
    wg = wg_task
    wg.name = 'test_wg_metadata_forwarded_to_submit'
    wg.metadata['label'] = 'Human-readable label'
    wg.submit(wait=True, timeout=30)
    assert wg.process.label == 'Human-readable label'


def test_wg_metadata_unset_falls_back_to_name(wg_task):
    """With no ``label`` in ``wg.metadata``, ``process.label`` falls back to the workgraph name
    exactly as it does today (pins current behavior)."""
    wg = wg_task
    wg.name = 'test_wg_metadata_unset_falls_back_to_name'
    assert 'label' not in wg.metadata
    wg.run()
    assert wg.process.label == 'test_wg_metadata_unset_falls_back_to_name'


def test_launch_metadata_overrides_wg_metadata(wg_task):
    """An explicit ``metadata={'label': ...}`` at launch time wins over ``wg.metadata['label']``."""
    wg = wg_task
    wg.name = 'test_launch_metadata_overrides_wg_metadata'
    wg.metadata['label'] = 'graph-level label'
    wg.run(metadata={'label': 'launch-time label'})
    assert wg.process.label == 'launch-time label'


def test_wg_metadata_merge_semantics(wg_task):
    """Launch-time metadata overrides matching keys; other ``wg.metadata`` keys survive."""
    wg = wg_task
    wg.name = 'test_wg_metadata_merge_semantics'
    wg.metadata['label'] = 'graph-level label'
    wg.metadata['description'] = 'graph-level description'
    wg.run(metadata={'label': 'launch-time label'})
    assert wg.process.label == 'launch-time label'
    assert wg.process.description == 'graph-level description'


@pytest.mark.parametrize(
    'metadata, offender',
    [
        ({'bad_key': 1}, 'bad_key'),  # a key outside the schema
        ({'store_provenance': 'no'}, 'store_provenance'),  # a declared key, wrong type under strict mode
    ],
)
def test_launch_metadata_override_validated(wg_task, metadata, offender):
    """A launch-time ``metadata`` override goes through the same schema as ``wg.metadata``."""
    wg = wg_task
    wg.name = 'test_launch_metadata_override_validated'
    with pytest.raises(ValueError, match=offender):
        wg.run(metadata=metadata)
    assert wg.process is None


def test_wg_metadata_does_not_affect_name_or_process_label(wg_task):
    """Setting ``wg.metadata`` leaves identity (``name``/``process_label``) untouched."""
    wg = wg_task
    wg.name = 'test_wg_metadata_does_not_affect_name_or_process_label'
    wg.metadata['label'] = 'some display label'
    wg.run()
    assert wg.name == 'test_wg_metadata_does_not_affect_name_or_process_label'
    assert wg.process.process_label == 'WorkGraph<test_wg_metadata_does_not_affect_name_or_process_label>'


def test_wg_metadata_roundtrips_through_dict(decorated_add):
    """``wg.metadata`` survives a ``to_dict``/``from_dict`` round trip, in the existing ``metadata`` slot."""
    wg = WorkGraph('test_wg_metadata_roundtrips_through_dict', metadata={'label': 'round-trip label'})
    wg.add_task(decorated_add, x=2, y=3)
    wgdata = wg.to_dict()
    assert wgdata['metadata']['label'] == 'round-trip label'
    wg2 = WorkGraph.from_dict(wgdata)
    assert wg2.metadata['label'] == 'round-trip label'
    assert wg2.name == 'test_wg_metadata_roundtrips_through_dict'


@pytest.mark.parametrize(
    'metadata, offender',
    [
        ({'bad_key': 1}, 'bad_key'),  # a key outside the schema
        ({'label': [1, 2, 3]}, 'label'),  # a declared key holding the wrong type
    ],
)
def test_wg_metadata_refused_on_construction(metadata, offender):
    """The constructor's ``metadata`` kwarg is validated straight away, naming the offender."""
    with pytest.raises(ValueError, match=offender):
        WorkGraph('test_wg_metadata_refused_on_construction', metadata=metadata)


@pytest.mark.parametrize(
    'field, bad_value',
    [
        ('store_provenance', 'no'),  # bool: lax mode would coerce each of these
        ('store_provenance', 0),
        ('store_provenance', 'yes'),
        ('pk', '12'),  # Optional[int]: a numeral spelled as a string
        ('label', 1),  # str: an int lax mode would stringify
    ],
)
def test_wg_metadata_strict_mode_refuses(field, bad_value):
    """Strict mode refuses values lax mode would have coerced into the field's type."""
    with pytest.raises(ValueError, match=field):
        WorkGraph.validate_metadata({field: bad_value})


@pytest.mark.parametrize('field, value', [('store_provenance', True), ('store_provenance', False), ('pk', 12)])
def test_wg_metadata_strict_mode_accepts(field, value):
    """A value of the field's own type still passes under strict mode, unchanged."""
    assert WorkGraph.validate_metadata({field: value})[field] == value


@pytest.fixture(params=['item', 'reassign'])
def wg_with_bad_key(request, wg_task):
    """``wg_task`` carrying a key outside the schema, set by item assignment or by replacing the dict."""
    wg = wg_task
    wg.name = request.node.name
    if request.param == 'item':
        wg.metadata['bad_key'] = 1
    else:
        wg.metadata = {'bad_key': 1}
    return wg


def test_wg_metadata_bad_key_refused_at_the_boundaries_not_the_line(wg_with_bad_key):
    """``wg.metadata`` is a plain dict: a bad key sticks at the line and is refused at every boundary.

    This is the design, not an oversight — validation lives at construction, load,
    serialization and launch, so nothing has to wrap the dict.
    """
    wg = wg_with_bad_key
    assert wg.metadata['bad_key'] == 1  # no complaint at the assignment
    with pytest.raises(ValueError, match='bad_key'):
        wg.to_dict()
    with pytest.raises(ValueError, match='bad_key'):
        wg.to_engine_inputs()
    with pytest.raises(ValueError, match='bad_key'):
        wg.run()
    assert wg.process is None
    wg.metadata = {'label': 'a valid replacement'}
    assert wg.to_dict()['metadata']['label'] == 'a valid replacement'


def test_engine_launch_keys_come_from_the_engine_spec():
    """The launch half of the schema is introspected from ``WorkGraphEngine``, not hand-listed.

    If AiiDA's own metadata port grows or loses a name, the schema follows it.
    """
    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    assert ENGINE_LAUNCH_KEYS == frozenset(WorkGraphEngine.spec().inputs['metadata'].keys())
    assert ENGINE_LAUNCH_KEYS <= frozenset(WorkGraphMetadata.__annotations__)


def test_wg_metadata_declared_keys_disjoint_from_bookkeeping():
    """AiiDA's launch-metadata port names never collide with node-graph's bookkeeping keys.

    This is the tripwire for the #812 collision: were a bookkeeping key ever to
    share a name with an AiiDA metadata port, ``to_engine_inputs`` could no
    longer tell which family a key belongs to, and a bookkeeping value could
    leak into the launch inputs. Checked against the two raw sources, so a
    collision cannot hide behind a set difference.
    """
    assert BOOKKEEPING_KEYS.isdisjoint(ENGINE_LAUNCH_KEYS)
    # sanity: the bookkeeping side is non-empty and known, not an accidental empty set
    assert BOOKKEEPING_KEYS == {'graph_class', 'definition', 'pk'}
    assert BOOKKEEPING_KEYS <= frozenset(WorkGraphMetadata.__annotations__)


def test_wg_metadata_task_graph_build_definition_key_not_rejected():
    """A nested ``@task.graph`` build writes node-graph's own ``definition`` bookkeeping key
    into the built graph's metadata (via ``graph.metadata.setdefault('definition', ...)``,
    node_graph/utils/graph.py). That write must not be rejected by `WorkGraph`'s validating
    metadata just because the key isn't an AiiDA launch key — it's declared bookkeeping too.
    """

    @task()
    def add(x, y):
        return x + y

    @task.graph()
    def inner_graph(x, y):
        return add(x, y).result

    with WorkGraph('test_wg_metadata_task_graph_build_definition_key_not_rejected') as wg:
        inner_graph(1, 2)
        wg.run()

    assert wg.process.is_finished_ok


def test_wg_metadata_legacy_wgdata_loads(decorated_add):
    """A ``wgdata`` whose ``metadata`` carries only bookkeeping keys (no AiiDA launch keys) —
    the shape written before this feature existed — still loads without error."""
    wg = WorkGraph('test_wg_metadata_legacy_wgdata_loads')
    wg.add_task(decorated_add, x=2, y=3)
    wgdata = wg.to_dict()
    assert set(wgdata['metadata']) <= {'graph_class', 'pk'}
    wg2 = WorkGraph.from_dict(wgdata)
    assert 'label' not in wg2.metadata
    assert wg2.name == 'test_wg_metadata_legacy_wgdata_loads'


def test_wg_metadata_unrecognized_legacy_key_raises_on_load(decorated_add):
    """A ``wgdata`` carrying a metadata key this schema has never declared — e.g.
    ``platform``/``worker_name``, written by an earlier aiida-workgraph version and
    by no current code — raises on load, naming the key and the graph, exactly as it
    would on fresh construction. A graph serialized with stray metadata keys must
    have them removed before it loads again.
    """
    wg = WorkGraph('test_wg_metadata_unrecognized_legacy_key_raises_on_load')
    wg.add_task(decorated_add, x=2, y=3)
    wgdata = wg.to_dict()
    wgdata['metadata']['worker_name'] = 'localhost'
    with pytest.raises(ValueError) as excinfo:
        WorkGraph.from_dict(wgdata)
    assert "graph 'test_wg_metadata_unrecognized_legacy_key_raises_on_load'" in str(excinfo.value)
    assert 'worker_name' in str(excinfo.value)


def test_wg_metadata_unset_graph_serializes_like_before(decorated_add):
    """A graph that never touches ``wg.metadata`` serializes exactly as it did before this feature."""
    wg = WorkGraph('test_wg_metadata_unset_graph_serializes_like_before')
    wg.add_task(decorated_add, x=2, y=3)
    wgdata = wg.to_dict()
    assert wgdata['metadata'] == {
        'graph_class': {'callable_name': 'WorkGraph', 'module_path': 'aiida_workgraph.workgraph'},
        'pk': None,
    }


def test_load_failure(create_process_node):
    node = create_process_node()
    with pytest.raises(ValueError, match=f'Process {node.pk} is not a WorkGraph'):
        WorkGraph.load(node.pk)


def test_organize_nested_inputs():
    """Merge sub properties to the root properties."""
    from .utils.test_workchain import WorkChainWithNestNamespace

    wg = WorkGraph('test_organize_nested_inputs')
    task1 = wg.add_task(WorkChainWithNestNamespace, name='task1')
    task1.set_inputs(
        {
            'add': {'x': '1'},
            'add.metadata': {
                'call_link_label': 'nest',
                'options': {'resources': {'num_cpus': 1}},
            },
            'add.metadata.options': {'resources': {'num_machines': 1}},
        }
    )
    inputs = wg.to_engine_inputs()
    data = {
        'metadata': {
            'call_link_label': 'nest',
            'options': {'resources': {'num_machines': 1}},
        },
        'x': '1',
    }
    assert inputs['tasks']['task1']['add'] == data


@pytest.mark.usefixtures('started_daemon_client')
def test_reset_message(wg_calcjob):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    from aiida.cmdline.utils.common import get_workchain_report

    wg = wg_calcjob
    wg.submit()
    timeout = 30
    wg.wait(tasks={'add1': ['RUNNING']}, timeout=timeout, interval=1)
    wg = WorkGraph.load(wg.process.pk)
    wg.tasks.add1.set_inputs({'y': orm.Int(10).store()})
    wg.save()
    wg.wait(timeout=timeout * 2)
    report = get_workchain_report(wg.process, 'REPORT')
    assert "Action: RESET. Tasks: ['add1']" in report


def test_restart_and_reset(wg_task):
    """Restart from a finished workgraph.
    Load the workgraph, modify the task, and restart the workgraph.
    Only the modified node and its child tasks will be rerun."""
    wg = wg_task
    wg.outputs.diff = wg.tasks.sumdiff1.outputs.diff
    wg.outputs.sum = wg.tasks.sumdiff2.outputs.sum
    wg.add_task(
        'workgraph.test_sum_diff',
        'sumdiff3',
        x=4,
        y=wg.tasks.sumdiff2.outputs.sum,
    )
    wg.name = 'test_restart_0'
    wg.run()
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.restart()
    wg1.name = 'test_restart_1'
    wg1.tasks.sumdiff2.set_inputs({'x': orm.Int(10).store()})
    wg1.run()
    assert wg1.tasks.sumdiff1.pk == wg.tasks.sumdiff1.pk
    assert wg1.tasks.sumdiff2.pk != wg.tasks.sumdiff2.pk
    assert wg1.tasks.sumdiff3.pk != wg.tasks.sumdiff3.pk
    assert wg1.tasks.sumdiff3.outputs.sum.value == 19
    wg1.reset()
    assert wg1.process is None
    assert wg1.tasks.sumdiff3.process is None
    assert wg1.tasks.sumdiff3.state == 'PLANNED'


@pytest.mark.skip(reason='This is break, opened ')
def test_extend_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph('test_graph_build')
    add1 = wg.add_task('workgraph.test_add', 'add1', x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group.build(x=0, y=4, z=5)
    # test wait
    add_multiply_wg.tasks.multiply.waiting_on.add('add')
    # extend workgraph
    wg.extend(add_multiply_wg, prefix='group_')
    assert 'group_add' in [task.name for task in wg.tasks.group_multiply.waiting_on]
    wg.add_link(add1.outputs[0], wg.tasks.group_add.inputs.x)
    wg.run()
    assert wg.tasks.group_multiply.outputs.result.value == 45


def test_workgraph_outputs(decorated_add):
    wg = WorkGraph('test_workgraph_outputs')
    wg.add_task(decorated_add, 'add1', x=2, y=3)
    wg.outputs.sum = wg.tasks.add1.outputs.result
    wg.run()
    assert wg.process.outputs.sum.value == 5
    # assert wg.process.outputs.add1.result.value == 5


@pytest.mark.usefixtures('started_daemon_client')
def test_wait_timeout(create_workgraph_process_node):
    wg = WorkGraph()
    wg.process = create_workgraph_process_node(state='running')
    with pytest.raises(
        TimeoutError,
        match='Timeout reached after 1 seconds while waiting for the WorkGraph:',
    ):
        wg.wait(timeout=1, interval=1)


def test_inputs_outputs(decorated_namespace_sum_diff):
    """Test the group inputs and outputs of the WorkGraph."""

    wg = WorkGraph(
        name='test_inputs_outputs',
        inputs=spec.namespace(x=Any, nested=spec.namespace(x=Any)),
    )
    wg.inputs = {'x': 1, 'nested.x': 2}
    # same as
    # wg.add_input("workgraph.any", "x")
    # wg.add_input("workgraph.namespace", "nested")
    # wg.add_input("workgraph.any", "nested.x")
    # wg.inputs.x = 1
    # wg.inputs.nested.x = 2
    wg.add_task(decorated_namespace_sum_diff, name='sum_diff1', x=wg.inputs.x, y=3)
    wg.tasks.sum_diff1.inputs.nested.x = wg.inputs.nested.x
    wg.tasks.sum_diff1.inputs.nested.y = 3
    wg.outputs.sum = wg.tasks.sum_diff1.outputs.sum
    wg.outputs.nested = {}
    wg.outputs.nested.sum = wg.tasks.sum_diff1.outputs.nested.sum
    # same as
    # wg.add_output("workgraph.namespace", "nested")
    # wg.add_output("workgraph.any", "nested.sum")
    wg.run()
    assert wg.outputs.sum.value == 4
    assert wg.outputs.nested.sum.value == 5


def test_inputs_run_submit_api():
    """Test running a WorkGraph with inputs provided in the `run` and `submit` APIs."""

    def generate_workgraph():
        with WorkGraph(inputs=spec.namespace(x=Any, y=Any)) as wg:
            wg.outputs.sum = wg.inputs.x + wg.inputs.y
        return wg

    wg = generate_workgraph()
    wg.run(inputs={'x': 1, 'y': 2})

    assert wg.outputs.sum.value == 3

    wg = generate_workgraph()
    wg.submit(inputs={'x': 3, 'y': 4}, wait=True)

    assert wg.outputs.sum.value == 7


def test_run_workgraph_builder():
    """Test running a WorkGraph using the WorkGraphEngine builder."""
    from aiida_workgraph.engine.workgraph import WorkGraphEngine
    from aiida.engine import run_get_node

    @task
    def add(x, y):
        """A simple task that adds two numbers."""
        return x + y

    wg = WorkGraph()
    wg.add_task(add, x=1, y=2)
    wgdata = wg.to_engine_inputs()
    builder = WorkGraphEngine.get_builder()
    builder._update(wgdata)
    _, node = run_get_node(builder)
    wg.process = node
    wg.update()
    assert wg.tasks.add.outputs.result.value == 3


def test_calling_workgraph_in_context_manager():
    """Test calling a `WorkGraph` in a context manager."""

    @task
    def add(x, y):
        return x + y

    with WorkGraph(inputs=spec.namespace(x=Any, y=Any), outputs=spec.namespace(sum=Any)) as wg1:
        add_outputs = add(x=wg1.inputs.x, y=wg1.inputs.y)  # add
        add1_outputs = add(x=add_outputs.result, y=1)
        wg1.outputs.sum = add1_outputs.result

    with WorkGraph() as wg2:
        sub_outputs = wg1({'x': 1, 'y': 2})
        add_outputs = add(x=sub_outputs.sum, y=5)
        wg2.outputs.sum = add_outputs.result

    wg2.run()

    assert wg2.outputs.sum.value == 9


def test_expose_task_spec():
    from aiida_workgraph import task
    from aiida_workgraph.socket_spec import namespace as ns

    @task()
    def test_calc(x: int) -> ns(square=int, double=int):
        return {'square': x * x, 'double': x + x}

    @task()
    def add_multiply(data: ns(x=int, y=int)) -> ns(sum=int, product=int):
        return {'sum': data['x'] + data['y'], 'product': data['x'] * data['y']}

    out = ns(out1=add_multiply.outputs, out2=test_calc.outputs['square'])

    @task.graph()
    def test_graph(x: int, data: ns(y=int)) -> out:
        am = add_multiply(data={'x': x, 'y': data['y']})
        tc = test_calc(x)
        return {'out1': am, 'out2': tc.square}

    wg = test_graph.build(x=1, data={'y': 2})
    wg.run()
    assert wg.outputs.out1.sum.value == 3
    assert wg.outputs.out1.product.value == 2
    assert wg.outputs.out2.value == 1


def test_update_outputs():
    """Test the update method of the WorkGraph."""
    from aiida_workgraph.orm.workgraph import WorkChainNode
    from aiida.common.links import LinkType

    def create_process_node(state='finished', exit_status=0, outputs: dict = None):
        node = WorkChainNode()
        node.store()
        # add outputs
        for key, value in outputs.items():
            if isinstance(value, dict):
                for k, v in value.items():
                    v.base.links.add_incoming(node, link_type=LinkType.RETURN, link_label=f'{key}__{k}')
            else:
                value.base.links.add_incoming(node, link_type=LinkType.RETURN, link_label=key)
        node.set_process_state(state)
        node.set_exit_status(exit_status)
        node.seal()
        return node

    # namespace output
    wg = WorkGraph('test_update_outputs', outputs=spec.namespace(y=Any))
    wg.process = create_process_node(outputs={'y': orm.Int(1).store()})
    wg.update()
    assert wg.state == 'FINISHED'
    assert wg.outputs.y.value == 1
    # dynamic output
    wg = WorkGraph('test_update_outputs', outputs=spec.dynamic(Any))
    wg.process = create_process_node(
        outputs={'x': orm.Int(0).store(), 'y': orm.Int(1).store(), 'z': orm.Int(2).store()}
    )
    wg.update()
    assert wg.outputs.x.value == 0
    assert wg.outputs.y.value == 1
    assert wg.outputs.z.value == 2
    # nested namespace output
    wg = WorkGraph('test_update_outputs', outputs=spec.namespace(x=Any, nested=spec.namespace(y=Any, z=Any)))
    wg.process = create_process_node(
        outputs={'x': orm.Int(0).store(), 'nested': {'y': orm.Int(1).store(), 'z': orm.Int(2).store()}}
    )
    wg.update()
    assert wg.outputs.x.value == 0
    assert wg.outputs.nested.y.value == 1
    assert wg.outputs.nested.z.value == 2
