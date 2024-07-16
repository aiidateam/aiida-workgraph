import pytest
from aiida_workgraph import WorkGraph
import os


@pytest.mark.usefixtures("started_daemon_client")
def test_submit(wg_calcjob: WorkGraph) -> None:
    """Submit simple calcjob."""
    wg = wg_calcjob
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True)
    os.system("verdi process list -a")
    os.system(f"verdi process report {wg.pk}")
    assert wg.tasks["add2"].outputs["sum"].value == 9
