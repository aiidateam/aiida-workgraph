import pytest
from aiida_workgraph import WorkGraph


def test_run(wg_calcfunction: WorkGraph) -> None:
    """Run simple calcfunction."""
    wg = wg_calcfunction
    wg.name = "test_run_calcfunction"
    wg.run()
    print("state: ", wg.state)
    # print("results: ", results[])
    assert wg.tasks["sumdiff2"].node.outputs.sum == 9
    assert wg.tasks["sumdiff2"].outputs["sum"].value == 9


@pytest.mark.usefixtures("started_daemon_client")
def test_submit(wg_calcfunction: WorkGraph) -> None:
    """Submit simple calcfunction."""
    wg = wg_calcfunction
    wg.name = "test_submit_calcfunction"
    wg.submit(wait=True)
    # print("results: ", results[])
    assert wg.tasks["sumdiff2"].outputs["sum"].value == 9
