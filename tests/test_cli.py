from click.testing import CliRunner
from aiida_workgraph.cli.cmd_workgraph import workgraph
import pytest
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida_workgraph import task

AddTask = task(ArithmeticAddCalculation)


@task.graph
def add_graph(x, y, code):
    return AddTask(code=code, x=x, y=y, metadata={'options': {'sleep': 15}}).sum


@task.graph
def two_add_graph(x, y, code):
    sum1 = AddTask(code=code, x=x, y=y, metadata={'options': {'sleep': 15}}).sum
    sum2 = AddTask(code=code, x=sum1, y=y, metadata={'options': {'sleep': 15}}).sum
    return sum2


def test_workgraph():
    """Test ``verdi group path ls``"""
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ['--help'])
    assert result.exit_code == 0, result.exception
    print(result.output)


def test_task():
    """Test ``verdi group path ls``"""
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ['task', '--help'])
    assert result.exit_code == 0, result.exception
    print(result.output)


@pytest.mark.usefixtures('started_daemon_client')
def test_task_pause_play(add_code):
    import time

    wg = add_graph.submit(1, 2, add_code)
    cli_runner = CliRunner()
    # wait the first task to be created/launched
    name = wg.tasks[-1].name
    print(f'Waiting for task {name} to be CREATED/WAITING')
    wg.wait(tasks={name: ['CREATED', 'RUNNING', 'WAITING']}, timeout=20, interval=2)
    result = cli_runner.invoke(workgraph, ['task', 'pause', str(wg.pk), name])
    assert result.exit_code == 0, result.exception
    tstart = time.time()
    node = wg.tasks[name].node
    # wait for the task to be paused
    print('state', node.process_state, node.process_status, node.paused)
    while not node.paused:
        print(node.process_state, node.process_status)
        time.sleep(1)
        if time.time() - tstart > 20:
            raise TimeoutError(f'Task {name} was not paused within 20s')
    # play the second task
    result = cli_runner.invoke(workgraph, ['task', 'play', str(wg.pk), name])
    assert result.exit_code == 0, result.exception
    print('state', node.process_state, node.process_status, node.paused)
    while node.paused:
        print(node.process_state, node.process_status)
        time.sleep(1)
        if time.time() - tstart > 20:
            raise TimeoutError(f'Task {name} was not paused within 20s')


@pytest.mark.usefixtures('started_daemon_client')
def test_task_kill(add_code, capsys):
    import time

    wg = add_graph.submit(1, 2, add_code)
    cli_runner = CliRunner()
    # wait the first task to be created/launched
    name = wg.tasks[-1].name
    print(f'Waiting for task {name} to be CREATED/WAITING')
    wg.wait(tasks={name: ['CREATED', 'RUNNING', 'WAITING']}, timeout=20, interval=2)
    result = cli_runner.invoke(workgraph, ['task', 'kill', str(wg.pk), name])
    assert result.exit_code == 0, result.exception
    tstart = time.time()
    node = wg.tasks[name].node
    while not node.is_killed:
        print(node.process_state, node.process_status)
        time.sleep(1)
        if time.time() - tstart > 20:
            raise TimeoutError(f'Task {name} was not killed within 20s')
    assert node.is_killed


@pytest.mark.usefixtures('started_daemon_client')
def test_task_list(add_code, capsys):
    wg = add_graph.build(1, 2, add_code)
    wg.save()
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ['task', 'list', str(wg.pk)])
    assert result.exit_code == 0, result.exception
    assert 'ArithmeticAddCalculation        PLANNED' in result.output
