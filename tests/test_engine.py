import logging
import time
import pytest
from aiida_workgraph import WorkGraph
from aiida.cmdline.utils.common import get_workchain_report


@pytest.mark.usefixtures('started_daemon_client')
def test_run_order(decorated_add) -> None:
    """Test the order.
    Tasks should run in parallel and only depend on the input tasks."""
    wg = WorkGraph(name='test_run_order')
    wg.add_task(decorated_add, 'add0', x=2, y=0)
    wg.add_task(decorated_add, 'add1', x=2, y=1)
    wg.submit(wait=True)
    report = get_workchain_report(wg.process, 'REPORT')
    assert 'tasks ready to run: add0,add1' in report


@pytest.mark.skip(reason='The test is not stable.')
def test_reset_node(wg_engine: WorkGraph) -> None:
    """Modify a node during the excution of a WorkGraph."""
    from aiida import orm

    wg = wg_engine
    wg.name = 'test_reset'
    wg.submit()
    time.sleep(15)
    wg.tasks['add3'].set_inputs({'y': orm.Int(10).store()})
    wg.save()
    wg.wait()
    wg.update()
    assert wg.tasks.add5.outputs.sum == 21
    assert wg.process.base.extras.get('_workgraph_queue_index') == 1
    assert len(wg.process.base.extras.get('_workgraph_queue')) == 1


@pytest.mark.usefixtures('started_daemon_client')
def test_max_number_jobs(add_code) -> None:
    from aiida_workgraph import WorkGraph
    from aiida.orm import Int
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    wg = WorkGraph('test_max_number_jobs')
    N = 3
    # Create N nodes
    for i in range(N):
        wg.add_task(ArithmeticAddCalculation, name=f'add{i}', x=Int(1), y=Int(1), code=add_code)
    # Set the maximum number of running jobs inside the WorkGraph
    wg.max_number_jobs = 2
    wg.submit(wait=True, timeout=40)
    report = get_workchain_report(wg.process, 'REPORT')
    assert 'tasks ready to run: add2' in report
    wg.tasks.add2.outputs.sum.value == 2


class _RecordingEngine:
    """Stand-in for a running engine that records the pause/kill calls it receives.

    ``_schedule_rpc`` calls the callback straight away, so a wrong keyword name
    raises ``TypeError`` here just as it would on the process event loop.
    """

    logger = logging.getLogger(__name__)
    pid = 0

    def __init__(self) -> None:
        self.killed: tuple | None = None
        self.paused: str | None = None

    def _schedule_rpc(self, callback, *args, **kwargs):
        return callback(*args, **kwargs)

    def kill(self, msg_text: str | None = None, force_kill: bool = False) -> None:
        self.killed = (msg_text, force_kill)

    def pause(self, msg_text: str | None = None) -> None:
        self.paused = msg_text


def test_message_receive_reads_plumpy_message_text() -> None:
    """The pause and kill RPCs carry the text plumpy's ``MessageBuilder`` writes."""
    from plumpy.process_comms import MessageBuilder
    from aiida_workgraph.engine.workgraph import WorkGraphEngine

    engine = _RecordingEngine()

    WorkGraphEngine.message_receive(engine, None, MessageBuilder.kill(text='killed by test'))
    assert engine.killed == ('killed by test', False)

    WorkGraphEngine.message_receive(engine, None, MessageBuilder.kill(text='forced', force_kill=True))
    assert engine.killed == ('forced', True)

    WorkGraphEngine.message_receive(engine, None, MessageBuilder.pause(text='paused by test'))
    assert engine.paused == 'paused by test'
