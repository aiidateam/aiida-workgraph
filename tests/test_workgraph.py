import pytest
from aiida_workgraph import WorkGraph, task, spec
from aiida import orm
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from typing import Any


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
    assert len(wg.tasks) == 5
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
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    wg = wg_calcfunction
    wg.add_task(decorated_add, name="add1", x=2, y=3)
    metadata = {
        "options": {
            "resources": {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 2,
            },
        }
    }
    wg.add_task(ArithmeticAddCalculation, name="add2", x=4, metadata=metadata)
    wg.name = "test_save_load"
    wg.save()
    assert wg.process.process_state.value.upper() == "CREATED"
    assert wg.process.process_label == "WorkGraph<test_save_load>"
    assert wg.process.label == "test_save_load"
    wg2 = WorkGraph.load(wg.process.pk)
    assert len(wg.tasks) == len(wg2.tasks)
    # check the executor of the decorated task
    callable = wg2.tasks.add1.get_executor()
    assert callable == wg.tasks.add1.get_executor()
    assert wg.tasks.add2.inputs.metadata._value == wg2.tasks.add2.inputs.metadata._value
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

    wg = WorkGraph("test_organize_nested_inputs")
    task1 = wg.add_task(WorkChainWithNestNamespace, name="task1")
    task1.set_inputs(
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
    assert inputs["tasks"]["task1"]["add"] == data


@pytest.mark.usefixtures("started_daemon_client")
def test_reset_message(wg_calcjob, capsys):
    """Modify a node and save the workgraph.
    This will add a message to the workgraph_queue extra field."""
    from aiida.cmdline.utils.common import get_workchain_report

    wg = wg_calcjob
    wg.submit()
    timeout = 30
    wg.wait(tasks={"add1": ["RUNNING"]}, timeout=timeout, interval=1)
    wg = WorkGraph.load(wg.process.pk)
    wg.tasks.add1.set_inputs({"y": orm.Int(10).store()})
    wg.save()
    wg.wait(timeout=timeout * 2)
    report = get_workchain_report(wg.process, "REPORT")
    assert "Action: RESET. Tasks: ['add1']" in report


def test_restart_and_reset(wg_calcfunction):
    """Restart from a finished workgraph.
    Load the workgraph, modify the task, and restart the workgraph.
    Only the modified node and its child tasks will be rerun."""
    wg = wg_calcfunction
    wg.outputs.diff = wg.tasks.sumdiff1.outputs.diff
    wg.outputs.sum = wg.tasks.sumdiff2.outputs.sum
    wg.add_task(
        "workgraph.test_sum_diff",
        "sumdiff3",
        x=4,
        y=wg.tasks.sumdiff2.outputs.sum,
    )
    wg.name = "test_restart_0"
    wg.run()
    wg1 = WorkGraph.load(wg.process.pk)
    wg1.restart()
    wg1.name = "test_restart_1"
    wg1.tasks.sumdiff2.set_inputs({"x": orm.Int(10).store()})
    wg1.run()
    assert wg1.tasks.sumdiff1.pk == wg.tasks.sumdiff1.pk
    assert wg1.tasks.sumdiff2.pk != wg.tasks.sumdiff2.pk
    assert wg1.tasks.sumdiff3.pk != wg.tasks.sumdiff3.pk
    assert wg1.tasks.sumdiff3.outputs.sum.value == 19
    wg1.reset()
    assert wg1.process is None
    assert wg1.tasks.sumdiff3.process is None
    assert wg1.tasks.sumdiff3.state == "PLANNED"


@pytest.mark.skip(reason="This is break, opened ")
def test_extend_workgraph(decorated_add_multiply_group):
    from aiida_workgraph import WorkGraph

    wg = WorkGraph("test_graph_build")
    add1 = wg.add_task("workgraph.test_add", "add1", x=2, y=3)
    add_multiply_wg = decorated_add_multiply_group.build(x=0, y=4, z=5)
    # test wait
    add_multiply_wg.tasks.multiply.waiting_on.add("add")
    # extend workgraph
    wg.extend(add_multiply_wg, prefix="group_")
    assert "group_add" in [task.name for task in wg.tasks.group_multiply.waiting_on]
    wg.add_link(add1.outputs[0], wg.tasks.group_add.inputs.x)
    wg.run()
    assert wg.tasks.group_multiply.outputs.result.value == 45


def test_workgraph_outputs(decorated_add):
    wg = WorkGraph("test_workgraph_outputs")
    wg.add_task(decorated_add, "add1", x=2, y=3)
    wg.outputs.sum = wg.tasks.add1.outputs.result
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


def test_inputs_outputs(decorated_namespace_sum_diff):
    """Test the group inputs and outputs of the WorkGraph."""

    wg = WorkGraph(
        name="test_inputs_outputs",
        inputs=spec.namespace(x=Any, nested=spec.namespace(x=Any)),
    )
    wg.inputs = {"x": 1, "nested.x": 2}
    # same as
    # wg.add_input("workgraph.any", "x")
    # wg.add_input("workgraph.namespace", "nested")
    # wg.add_input("workgraph.any", "nested.x")
    # wg.inputs.x = 1
    # wg.inputs.nested.x = 2
    wg.add_task(decorated_namespace_sum_diff, name="sum_diff1", x=wg.inputs.x, y=3)
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
    wg.run(inputs={"x": 1, "y": 2})

    assert wg.outputs.sum.value == 3

    wg = generate_workgraph()
    wg.submit(inputs={"x": 3, "y": 4}, wait=True)

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
    wgdata = wg.prepare_inputs()
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

    with WorkGraph(inputs=spec.namespace(x=Any, y=Any)) as wg1:
        add_outputs = add(x=wg1.inputs.x, y=wg1.inputs.y)  # add
        add1_outputs = add(x=add_outputs.result, y=1)
        wg1.outputs.sum = add1_outputs.result

    with WorkGraph() as wg2:
        sub_outputs = wg1({"x": 1, "y": 2})
        add_outputs = add(x=sub_outputs.sum, y=5)
        wg2.outputs.sum = add_outputs.result

    wg2.run()

    assert wg2.outputs.sum.value == 9


def test_expose_task_spec():
    from aiida_workgraph import task
    from aiida_workgraph.socket_spec import namespace as ns

    @task()
    def test_calc(x: int) -> ns(square=int, double=int):
        return {"square": x * x, "double": x + x}

    @task()
    def add_multiply(data: ns(x=int, y=int)) -> ns(sum=int, product=int):
        return {"sum": data["x"] + data["y"], "product": data["x"] * data["y"]}

    out = ns(out1=add_multiply.outputs, out2=test_calc.outputs["square"])

    @task.graph()
    def test_graph(x: int, data: ns(y=int)) -> out:
        am = add_multiply(data={"x": x, "y": data["y"]})
        tc = test_calc(x)
        return {"out1": am, "out2": tc.square}

    wg = test_graph.build(x=1, data={"y": 2})
    wg.run()
    assert wg.outputs.out1.sum.value == 3
    assert wg.outputs.out1.product.value == 2
    assert wg.outputs.out2.value == 1
