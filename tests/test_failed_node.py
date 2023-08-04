import aiida

aiida.load_profile()


def test_failed_node(decorated_sqrt, decorated_add):
    """Submit simple calcfunction."""
    from aiida_worktree import WorkTree
    from aiida.orm import Float

    wt = WorkTree(name="test_failed_node")
    wt.nodes.new(decorated_add, "add1", x=Float(1), y=Float(2))
    wt.nodes.new(decorated_sqrt, "sqrt1", x=Float(-1))
    wt.submit(wait=True)
    # print("results: ", results[])
    assert wt.process.exit_status == 302
