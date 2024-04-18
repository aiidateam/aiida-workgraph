import aiida
import time
import pytest

aiida.load_profile()


def test_run_order(wt_engine):
    """Test the order.
    Nodes should run in parallel and only depend on the input nodes."""
    wt = wt_engine
    wt.submit(wait=True)
    wt.nodes["add2"].ctime < wt.nodes["add4"].ctime


@pytest.mark.skip(reason="The test is not stable.")
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


def test_max_number_jobs():
    from aiida_worktree import WorkTree, build_node
    from aiida.orm import load_code
    from aiida.orm import Int

    # Use the calcjob: ArithmeticAddCalculation
    arithmetic_add = build_node(
        "aiida.calculations.arithmetic.add.ArithmeticAddCalculation"
    )
    code = load_code("add@localhost")

    wt = WorkTree("test_max_number_jobs")
    N = 15
    # Create N nodes
    for i in range(N):
        temp = wt.nodes.new(
            arithmetic_add, name=f"add{i}", x=Int(1), y=Int(1), code=code
        )
        # Set a sleep option for each job (e.g., 2 seconds per job)
        temp.set({"metadata.options.sleep": 2})

    # Set the maximum number of running jobs inside the WorkTree
    wt.max_number_jobs = 5
    wt.submit(wait=True)
    wt.nodes["add2"].ctime < wt.nodes["add10"].ctime
