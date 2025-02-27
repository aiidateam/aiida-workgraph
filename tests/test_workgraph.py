import pytest
from aiida_workgraph import WorkGraph
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


def test_from_dict(decorated_add):
    """Export NodeGraph to dict."""
    wg = WorkGraph("test_from_dict")
    task1 = wg.add_task(decorated_add, x=2, y=3)
    wg.add_task("workgraph.test_sum_diff", name="sumdiff2", x=4, y=task1.outputs.result)
    wgdata = wg.to_dict()
    wg1 = WorkGraph.from_dict(wgdata)
    assert len(wg.tasks) == len(wg1.tasks)
    assert len(wg.links) == len(wg1.links)


def test_add_task():
    """Add add task."""
    wg = WorkGraph("test_add_task")
    add1 = wg.add_task(ArithmeticAddCalculation, name="add1")
    add2 = wg.add_task(ArithmeticAddCalculation, name="add2")
    wg.add_link(add1.outputs.sum, add2.inputs.x)
    assert len(wg.tasks) == 2
    assert len(wg.links) == 1


def test_show_state(wg_calcfunction):
    from io import StringIO
    import sys

    # Redirect stdout to capture prints
    captured_output = StringIO()
    sys.stdout = captured_output
    # Call the method
    wg_calcfunction.name = "test_show_state"
    wg_calcfunction.show()
    # Reset stdout
    sys.stdout = sys.__stdout__
    # Check the output
    output = captured_output.getvalue()
    assert "WorkGraph: test_show_state, PK: None, State: CREATED" in output
    assert "sumdiff1" in output
    assert "PLANNED" in output


def test_save_load(wg_calcfunction, decorated_add):
    """Save the workgraph"""
    from aiida_workgraph.orm.pickled_function import PickledFunction

    wg = wg_calcfunction
    wg.add_task(decorated_add, name="add1", x=2, y=3)
    wg.name = "test_save_load"
    wg.save()
    assert wg.process.process_state.value.upper() == "CREATED"
    assert wg.process.process_label == "WorkGraph<test_save_load>"
    assert wg.process.label == "test_save_load"
    wg2 = WorkGraph.load(wg.process.pk)
    assert len(wg.tasks) == len(wg2.tasks)
    # check the executor of the decorated task
    callable = wg2.tasks.add1.get_executor()["callable"]
    assert isinstance(callable, PickledFunction)
    assert callable.value(1, 2) == 3
    # TODO, the following code is not working
    # wg2.save()
    # assert wg2.tasks.add1.executor == decorated_add
    # remove the extra, will raise an error


def test_load_failure(create_process_node):
    node = create_process_node()
    with pytest.raises(ValueError, match=f"Process {node.pk} is not a WorkGraph"):
        WorkGraph.load(node.pk)


def test_organize_nested_inputs():
    """Merge sub properties to the root properties."""
    from .utils.test_workchain import WorkChainWithNestNamespace
    from node_graph.utils import collect_values_inside_namespace

    wg = WorkGraph("test_organize_nested_inputs")
    task1 = wg.add_task(WorkChainWithNestNamespace, name="task1")
    task1.set(
        {
            "add": {"x": "1"},
            "add.metadata": {
                "call_link_label": "nest",
                "options": {"resources": {"num_cpus": 1}},
            },
            "add.metadata.options": {"resources": {"num_machines": 1}},
        }
    )
    inputs = wg.prepare_inputs()
    data = {
        "metadata": {
            "call_link_label": "nest",
            "options": {"resources": {"num_machines": 1}},
        },
        "x": "1",
    }
    collected_data = collect_values_inside_namespace(
        inputs["workgraph_data"]["tasks"]["task1"]["inputs"]["add"]
    )
    assert collected_data == data


@pytest.mark.usefixtures("started_daemon_client")
def test_reset_message(wg_calcjob):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    from aiida.cmdline.utils.common import get_workchain_report

    wg = wg_calcjob
    wg.submit()
    wg.wait(tasks={"add1": ["RUNNING"]}, timeout=30, interval=1)
    wg = WorkGraph.load(wg.process.pk)
    wg.tasks.add1.set({"y": orm.Int(10).store()})
    wg.save()
    wg.wait(timeout=30)
    report = get_workchain_report(wg.process, "REPORT")
    assert "Action: reset. {'add1'}" in report


def test_restart_and_reset(wg_calcfunction):
    """Restart from a finished workgraph.
    Load the workgraph, modify the task, and restart the workgraph.
    Only the modified node and its child tasks will be rerun."""
    wg = wg_calcfunction
    wg.add_task(
        "workgraph.test_sum_diff",
        "sumdiff3",
        x=4,
        y=wg.tasks["sumdiff2"].outputs.sum,
    )
    wg.name = "test_restart_0"
    wg.run()
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.restart()
    wg1.name = "test_restart_1"
    wg1.tasks["sumdiff2"].set({"x": orm.Int(10).store()})
    wg1.run()
    assert wg1.tasks["sumdiff1"].node.pk == wg.tasks["sumdiff1"].pk
    assert wg1.tasks["sumdiff2"].node.pk != wg.tasks["sumdiff2"].pk
    assert wg1.tasks["sumdiff3"].node.pk != wg.tasks["sumdiff3"].pk
    assert wg1.tasks["sumdiff3"].node.outputs.sum == 19
    wg1.reset()
    assert wg1.process is None
    assert wg1.tasks["sumdiff3"].process is None
    assert wg1.tasks["sumdiff3"].state == "PLANNED"


def test_extend_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph("test_graph_build")
    add1 = wg.add_task("workgraph.test_add", "add1", x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group(x=0, y=4, z=5)
    # test wait
    add_multiply_wg.tasks["multiply1"].waiting_on.add("add1")
    # extend workgraph
    wg.extend(add_multiply_wg, prefix="group_")
    assert "group_add1" in [
        task.name for task in wg.tasks["group_multiply1"].waiting_on
    ]
    wg.add_link(add1.outputs[0], wg.tasks["group_add1"].inputs.x)
    wg.run()
    assert wg.tasks["group_multiply1"].node.outputs.result == 45


def test_workgraph_group_outputs(decorated_add):
    wg = WorkGraph("test_workgraph_group_outputs")
    wg.add_task(decorated_add, "add1", x=2, y=3)
    wg.group_outputs = [
        {"name": "sum", "from": "add1.result"},
        # {"name": "add1", "from": "add1"},
    ]
    wg.run()
    assert wg.process.outputs.sum.value == 5
    # assert wg.process.outputs.add1.result.value == 5


@pytest.mark.usefixtures("started_daemon_client")
def test_wait_timeout(create_workgraph_process_node):
    wg = WorkGraph()
    wg.process = create_workgraph_process_node(state="running")
    with pytest.raises(
        TimeoutError,
        match="Timeout reached after 1 seconds while waiting for the WorkGraph:",
    ):
        wg.wait(timeout=1, interval=1)

