import pytest
from aiida_workgraph import WorkGraph
from aiida_workgraph.tools.pytest_fixtures.scheduler import DEFAULT_TEST_SCHEDULER_NAME
from aiida_workgraph.engine.scheduler import Scheduler
from aiida_workgraph.orm.scheduler import SchedulerNode
from aiida.orm import CalcJobNode, WorkflowNode

from aiida_workgraph.engine.scheduler.client import get_scheduler_node
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
import datetime


def test_scheduler_node_init() -> None:
    """Test Scheduler node."""
    s = SchedulerNode(name="test", max_calcjobs=1, max_workflows=2, max_processes=3)
    assert s.name == "test"
    assert s.max_calcjobs == 1
    assert s.max_workflows == 2
    assert s.max_processes == 3


def test_scheduler_node_append_process() -> None:
    """Test Scheduler node."""
    s = SchedulerNode(name="test")
    calc_node = CalcJobNode().store()
    wf_node = WorkflowNode().store()
    s.append_running_process(calc_node.pk)
    s.append_running_process(wf_node.pk)
    assert s.running_process == [calc_node.pk, wf_node.pk]
    assert s.running_calcjob == [calc_node.pk]
    assert s.running_workflow == [wf_node.pk]
    # remove calcjob
    s.remove_running_process(calc_node.pk)
    assert s.running_process == [wf_node.pk]
    assert s.running_calcjob == []
    # remove workflow
    s.remove_running_process(wf_node.pk)
    assert s.running_process == []
    assert s.running_workflow == []


def test_scheduler_init() -> None:
    """Test Scheduler."""
    s = Scheduler(name="test", max_calcjobs=1, max_workflows=2, max_processes=3)
    assert s.name == "test"
    assert s.node.max_calcjobs == 1
    assert s.node.max_workflows == 2
    assert s.node.max_processes == 3


@pytest.mark.usefixtures("started_scheduler_client")
def test_scheduler_set_limit() -> None:
    from aiida_workgraph.engine.scheduler.client import get_scheduler_node

    s = get_scheduler_node(name=DEFAULT_TEST_SCHEDULER_NAME)
    response = Scheduler.set_max_calcjobs(DEFAULT_TEST_SCHEDULER_NAME, max_calcjobs=5)
    # wait for scheduler to be updated
    response.result().result()
    assert s.max_calcjobs == 5
    response = Scheduler.set_max_workflows(DEFAULT_TEST_SCHEDULER_NAME, max_workflows=6)
    response.result().result()
    assert s.max_workflows == 6
    response = Scheduler.set_max_processes(DEFAULT_TEST_SCHEDULER_NAME, max_processes=7)
    response.result().result()
    assert s.max_processes == 7


def test_scheduler_consume_process() -> None:
    """Test Scheduler consume process."""
    s = Scheduler(name="test2", max_calcjobs=1, max_workflows=1, max_processes=5)
    assert s.node.max_calcjobs == 1
    assert s.node.max_workflows == 1
    assert s.node.max_processes == 5
    calc_node = CalcJobNode().store()
    wf_node = WorkflowNode().store()
    assert s._has_capacity(calc_node) is True
    assert s._has_capacity(wf_node) is True
    s.node.append_running_process(calc_node.pk)
    s.node.append_running_process(wf_node.pk)
    assert len(s.node.running_calcjob) == 1
    assert len(s.node.running_workflow) == 1
    assert SchedulerNode.is_top_level_workflow(wf_node) is True
    assert s._has_capacity(calc_node) is False
    assert s._has_capacity(wf_node) is False


@pytest.mark.usefixtures("started_daemon_client")
@pytest.mark.usefixtures("started_scheduler_client")
def test_submit(add_code) -> None:
    """Submit simple calcjob."""
    wg = WorkGraph(name="test_debug_math")
    add1 = wg.add_task(ArithmeticAddCalculation, "add1", x=2, y=3, code=add_code)
    add2 = wg.add_task(ArithmeticAddCalculation, "add2", x=4, y=3, code=add_code)
    add1.set({"metadata.options.sleep": 5})
    add2.set({"metadata.options.sleep": 5})
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True, timeout=30, scheduler=DEFAULT_TEST_SCHEDULER_NAME)
    wg.update()
    assert add2.process.mtime - add1.process.mtime < datetime.timedelta(seconds=5)

    #
    status = Scheduler.get_status(name=DEFAULT_TEST_SCHEDULER_NAME)
    assert status is not None
    assert "running_process" in status
    status = Scheduler.set_max_calcjobs(
        name=DEFAULT_TEST_SCHEDULER_NAME, max_calcjobs=1
    )
    scheduler = get_scheduler_node(name=DEFAULT_TEST_SCHEDULER_NAME)
    assert scheduler.max_calcjobs == 1
    #
    wg.reset()
    wg.submit(wait=True, timeout=30, scheduler=DEFAULT_TEST_SCHEDULER_NAME)
    wg.update()
    assert add2.process.mtime - add1.process.mtime > datetime.timedelta(seconds=5)
