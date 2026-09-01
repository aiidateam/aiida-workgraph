import pytest
from aiida_workgraph import task
from typing import Annotated
from node_graph.task_spec import SchemaSource


def test_create_task_from_calcJob(add_code) -> None:
    """Test creating a task from a CalcJob."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    AddTask = task()(ArithmeticAddCalculation)
    metadata = {
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 2,
            },
        }
    }

    @task.graph
    def test_calcjob(inputs: Annotated[dict, AddTask.inputs]) -> Annotated[dict, AddTask.outputs]:
        return AddTask(x=inputs['x'], y=inputs['y'], code=inputs['code'], metadata=metadata)

    _, wg = test_calcjob.run_get_graph({'x': 2, 'y': 3, 'code': add_code, 'metadata': metadata})

    assert wg.outputs.sum.value == 5
    assert wg.tasks[-1].spec.schema_source == SchemaSource.CALLABLE
    assert wg.tasks[-1].get_executor().callable == ArithmeticAddCalculation


def test_calcjob_port_help_reaches_sockets() -> None:
    """The help text declared on AiiDA ports is carried onto the task's input sockets."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
    from aiida_workgraph import WorkGraph

    wg = WorkGraph()
    add1 = wg.add_task(ArithmeticAddCalculation)
    port_help = ArithmeticAddCalculation.spec().inputs['x'].help
    assert port_help
    assert add1.inputs.x._metadata.help == port_help


def test_missing_calcjob_inputs_carry_port_help() -> None:
    """Missing-input entries for a wrapped CalcJob expose the port help text."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
    from aiida_workgraph import WorkGraph
    from aiida_workgraph.errors import MissingRequiredInputsError

    wg = WorkGraph()
    wg.add_task(ArithmeticAddCalculation, name='add1')
    with pytest.raises(MissingRequiredInputsError) as excinfo:
        wg.run()
    entries = {entry.socket_path: entry for entry in excinfo.value.missing}
    assert {'add1.x', 'add1.y'} <= set(entries)
    assert entries['add1.x'].help == ArithmeticAddCalculation.spec().inputs['x'].help
