from aiida.cmdline.utils.common import get_workchain_report
import aiida
import time

aiida.load_profile()


def test_run_order(wt_engine):
    """Test the order.
    Nodes should run in parallel and only depend on the input nodes."""
    wt = wt_engine
    wt.submit(wait=True)
    # get the report
    report = get_workchain_report(wt.process, "REPORT")
    lines = report.splitlines()
    steps = {}
    for i, line in enumerate(lines):
        if "Run node:" in line:
            print(line)
            name = line.split("Run node:")[1].split(",")[0].strip()
            steps[name] = i
    # first run add0, add1
    # then add3
    # then add4
    # then add2
    # finally add5
    assert steps["add2"] > steps["add4"]


def test_reset_node(wt_engine):
    """Modify a node during the excution of a WorkTree."""
    wt = wt_engine
    wt.name = "test_reset"
    wt.submit()
    time.sleep(15)
    wt.nodes["add3"].set({"y": aiida.orm.Int(10).store()})
    wt.save()
    wt.wait()
    wt.update()
    assert wt.nodes["add5"].node.outputs.sum == 21
    assert wt.process.base.extras.get("worktree_queue_index") == 1
    assert len(wt.process.base.extras.get("worktree_queue")) == 1
