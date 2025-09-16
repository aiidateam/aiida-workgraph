import pytest
from aiida_workgraph import WorkGraph, task
from aiida_workgraph.socket import TaskSocketNamespace
from typing import Callable, Any, Annotated
from aiida_workgraph.manager import get_current_graph, set_current_graph
from aiida_workgraph import socket_spec as spec


def test_custom_outputs():
    """Test custom outputs."""

    @task
    def add_multiply(x, y) -> Annotated[dict, spec.namespace(sum=Any, product=Any)]:
        return {'sum': x + y, 'product': x * y}

    @task.graph(outputs=add_multiply.outputs)
    def test_graph(x, y):
        return add_multiply(x, y)

    wg = test_graph.build(1, 2)
    wg.run()
    assert wg.outputs.sum.value == 3
    assert wg.outputs.product.value == 2


def test_decorators_args() -> None:
    @task()
    def test(a, b=1, **c):
        print(a, b, c)

    n = test._spec.to_node()
    assert n.args_data['args'] == []
    assert set(n.args_data['kwargs']) == set(
        {
            'metadata',
            'function_data',
            'a',
            'process_label',
            'b',
            'function_inputs',
        }
    )
    assert n.args_data['var_args'] is None
    assert n.args_data['var_kwargs'] == 'c'
    assert set(n.get_output_names()) == {'_outputs', '_wait'}
    assert isinstance(n.inputs.c, TaskSocketNamespace)


def test_decorators_calcfunction_args() -> None:
    @task.calcfunction()
    def test(a, b=1, **c):
        print(a, b, c)

    metadata_kwargs = {f'{key}' for key in test._func.process_class.spec().inputs.ports['metadata'].ports.keys()}
    kwargs = set(test._func.process_class.spec().inputs.ports.keys())
    n = test._spec.to_node()
    assert n.args_data['args'] == []
    assert set(n.args_data['kwargs']) == set(kwargs)
    assert n.args_data['var_args'] is None
    assert n.args_data['var_kwargs'] == 'c'
    assert set(n.get_output_names()) == {'_outputs', '_wait'}
    assert isinstance(n.inputs.c, TaskSocketNamespace)
    assert set(n.inputs.metadata._get_keys()) == metadata_kwargs


@pytest.fixture(params=['decorator_factory', 'decorator'])
def task_function(request):
    if request.param == 'decorator_factory':

        @task()
        def test(a, b=1, **c):
            print(a, b, c)

    elif request.param == 'decorator':

        @task
        def test(a, b=1, **c):
            print(a, b, c)

    else:
        raise ValueError(f'{request.param} not supported.')
    return test


def test_decorators_task_args(task_function):
    n = task_function._spec.to_node()
    assert n.args_data['args'] == []
    assert set(n.args_data['kwargs']) == {
        'metadata',
        'function_data',
        'a',
        'process_label',
        'b',
        'function_inputs',
    }
    assert n.args_data['var_args'] is None
    assert n.args_data['var_kwargs'] == 'c'


@pytest.fixture(params=['decorator_factory', 'decorator'])
def task_workfunction(request):
    if request.param == 'decorator_factory':

        @task.workfunction()
        def test(a, b=1, **c):
            print(a, b, c)

    elif request.param == 'decorator':

        @task.workfunction
        def test(a, b=1, **c):
            print(a, b, c)

    else:
        raise ValueError(f'{request.param} not supported.')
    return test


def test_decorators_workfunction_args(task_workfunction) -> None:
    metadata_kwargs = {
        f'{key}' for key in task_workfunction._func.process_class.spec().inputs.ports['metadata'].ports.keys()
    }
    kwargs = set(task_workfunction._func.process_class.spec().inputs.ports.keys())
    #
    n = task_workfunction._spec.to_node()
    assert n.args_data['args'] == []
    assert set(n.args_data['kwargs']) == set(kwargs)
    assert n.args_data['var_args'] is None
    assert n.args_data['var_kwargs'] == 'c'
    assert set(n.get_output_names()) == {'_outputs', '_wait'}
    assert set(n.inputs.metadata._get_keys()) == metadata_kwargs


def test_decorators_parameters() -> None:
    """Test passing parameters to decorators."""

    @task.calcfunction(
        outputs=['sum', 'product'],
    )
    def test(a, b=1, **c):
        return {'sum': a + b, 'product': a * b}

    test1 = test._spec.to_node()
    assert test1.inputs['c']._link_limit == 1000000
    assert 'sum' in test1.get_output_names()
    assert 'product' in test1.get_output_names()


@pytest.fixture(params=['decorator_factory', 'decorator'])
def task_graph_task(request):
    if request.param == 'decorator_factory':

        @task.graph()
        def add_multiply_group(a, b=1, **c):
            pass

    elif request.param == 'decorator':

        @task.graph
        def add_multiply_group(a, b=1, **c):
            pass

    else:
        raise ValueError(f'{request.param} not supported.')

    return add_multiply_group


def test_decorators_graph_args(task_graph_task) -> None:
    # assert task_graph_task.identifier == "add_multiply_group"
    n = task_graph_task._spec.to_node()
    assert n.args_data['args'] == []
    assert n.args_data['kwargs'] == ['a', 'b']
    assert n.args_data['var_args'] is None
    assert n.args_data['var_kwargs'] == 'c'
    assert set(n.get_output_names()) == {'_outputs', '_wait'}


def test_inputs_outputs_workchain() -> None:
    from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

    wg = WorkGraph()
    task = wg.add_task(MultiplyAddWorkChain)
    assert 'metadata' in task.get_input_names()
    assert 'call_link_label' in task.inputs.metadata._get_keys()
    assert 'result' in task.get_output_names()


@pytest.mark.usefixtures('started_daemon_client')
def test_decorator_calcfunction(decorated_add: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name='test_decorator_calcfunction')
    wg.add_task(decorated_add, 'add1', x=2, y=3)
    wg.submit(wait=True, timeout=100)
    assert wg.tasks.add1.outputs.result.value == 5


def test_decorator_workfunction(decorated_add_multiply: Callable) -> None:
    """Run simple calcfunction."""

    wg = WorkGraph(name='test_decorator_workfunction')
    wg.add_task(decorated_add_multiply, 'add_multiply1', x=2, y=3, z=4)
    wg.run()
    assert wg.tasks['add_multiply1'].outputs.result.value == 20


def test_decorator_graph_namespace_outputs(decorated_add: Callable) -> None:
    """=Test namespace outputs in graph builder."""
    from aiida_workgraph.socket import TaskSocketNamespace, TaskSocket

    @task.graph
    def add_group(x, y, z) -> spec.namespace(add1=spec.dynamic(Any), sum=Any):
        outputs1 = decorated_add(x=x, y=y)
        return {'add1': outputs1, 'sum': outputs1.result}

    wg = WorkGraph()
    wg.add_task(add_group, x=2, y=3, z=4)
    assert isinstance(wg.tasks.add_group.outputs.add1, TaskSocketNamespace)
    assert isinstance(wg.tasks.add_group.outputs.sum, TaskSocket)


@pytest.mark.usefixtures('started_daemon_client')
def test_decorator_graph(decorated_add_multiply_group: Callable) -> None:
    """Test graph build."""
    wg = WorkGraph('test_graph')
    add1 = wg.add_task('workgraph.test_add', 'add1', x=2, y=3)
    add_multiply1 = wg.add_task(decorated_add_multiply_group, 'add_multiply1', y=3, z=4)
    sum_diff1 = wg.add_task('workgraph.test_sum_diff', 'sum_diff1')
    wg.add_link(add1.outputs[0], add_multiply1.inputs.x)
    wg.add_link(add_multiply1.outputs.result, sum_diff1.inputs.x)
    # use run to check if graph builder workgraph can be submit inside the engine
    wg.run()
    assert wg.tasks['add_multiply1'].process.outputs.result.value == 32
    assert wg.tasks['add_multiply1'].outputs.result.value == 32
    assert wg.tasks['sum_diff1'].outputs.sum.value == 32


def test_get_current_graph():
    g = get_current_graph()
    assert isinstance(g, WorkGraph)


def test_set_current_graph():
    @task()
    def add(x, y):
        return x + y

    sum = add(1, 2)
    g = get_current_graph()
    assert g == sum._node.graph
    g2 = WorkGraph()
    set_current_graph(g2)
    assert get_current_graph() == g2
