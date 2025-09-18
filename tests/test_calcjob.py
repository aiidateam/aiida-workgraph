from aiida_workgraph import task
from typing import Annotated


def test_create_task_from_calcJob(add_code) -> None:
    """Test creating a task from a CalcJob."""
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    AddTask = task()(ArithmeticAddCalculation)
    metadata = {
        "options": {
            "resources": {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 2,
            },
        }
    }

    @task.graph
    def test_calcjob(
        inputs: Annotated[dict, AddTask.inputs]
    ) -> Annotated[dict, AddTask.outputs]:
        return AddTask(
            x=inputs["x"], y=inputs["y"], code=inputs["code"], metadata=metadata
        )

    _, wg = test_calcjob.run_get_graph(
        {"x": 2, "y": 3, "code": add_code, "metadata": metadata}
    )

    assert wg.outputs.sum.value == 5
    assert wg.tasks[-1]._spec.mode == "decorator_build"
    assert wg.tasks[-1].get_executor().callable == ArithmeticAddCalculation
