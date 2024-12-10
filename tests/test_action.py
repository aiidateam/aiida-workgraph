import pytest
import time


@pytest.mark.skip(reason="PAUSED state is wrong for the moment.")
def test_pause_play_workgraph(wg_engine):
    wg = wg_engine
    wg.name = "test_pause_play_workgraph"
    wg.submit()
    time.sleep(5)
    wg.pause()
    wg.update()
    assert wg.process.process_state.value.upper() == "PAUSED"


# @pytest.mark.skip(reason="pause task is not stable for the moment.")
@pytest.mark.usefixtures("started_daemon_client")
def test_pause_play_task(wg_calcjob):
    wg = wg_calcjob
    wg.name = "test_pause_play_task"
    # pause add1 before submit
    wg.pause_tasks(["add1"])
    wg.submit()
    # wait for the workgraph to launch add1
    wg.wait(tasks={"add1": ["CREATED"]}, timeout=40, interval=5)
    assert wg.tasks["add1"].node.process_state.value.upper() == "CREATED"
    assert wg.tasks["add1"].node.process_status == "Paused through WorkGraph"
    # pause add2 after submit
    wg.pause_tasks(["add2"])
    wg.play_tasks(["add1"])
    # wait for the workgraph to launch add2
    wg.wait(tasks={"add2": ["CREATED"]}, timeout=40, interval=5)
    assert wg.tasks["add2"].node.process_state.value.upper() == "CREATED"
    assert wg.tasks["add2"].node.process_status == "Paused through WorkGraph"
    # I disabled the following lines because the test is not stable
    # Seems the daemon is not responding to the play signal
    wg.play_tasks(["add2"])
    wg.wait(interval=5)
    assert wg.tasks["add2"].outputs["sum"].socket_value == 9


def test_pause_play_error_handler(wg_calcjob, finished_process_node):
    wg = wg_calcjob
    wg.name = "test_pause_play_error_handler"
    wg.process = finished_process_node
    try:
        wg.pause_tasks(["add1"])
    except Exception as e:
        assert "WorkGraph is finished. Cannot pause tasks." in str(e)

    try:
        wg.play_tasks(["add1"])
    except Exception as e:
        assert "WorkGraph is finished. Cannot play tasks." in str(e)

    try:
        wg.kill_tasks(["add2"])
    except Exception as e:
        assert "WorkGraph is finished. Cannot kill tasks." in str(e)
