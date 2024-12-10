import pytest
from aiida_workgraph import WorkGraph


@pytest.mark.usefixtures("started_daemon_client")
def test_submit(wg_calcjob: WorkGraph) -> None:
    """Submit simple calcjob."""
    wg = wg_calcjob
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True)
    assert wg.tasks["add2"].outputs["sum"].socket_value == 9
