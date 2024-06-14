import aiida
from aiida_workgraph import WorkGraph


aiida.load_profile()


def test_submit(wg_calcjob: WorkGraph) -> None:
    """Submit simple calcjob."""
    wg = wg_calcjob
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True)
    assert wg.tasks["add2"].outputs["sum"].value == 9
