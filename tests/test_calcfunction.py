import aiida

aiida.load_profile()


def test_run(nt_calcfunction):
    """Run simple calcfunction."""
    nt = nt_calcfunction
    nt.name = "test_run_calcfunction"
    nt.run()
    print("state: ", nt.state)
    # print("results: ", results[])
    assert nt.nodes["sumdiff2"].node.outputs.sum == 9


def test_submit(nt_calcfunction):
    """Submit simple calcfunction."""
    nt = nt_calcfunction
    nt.name = "test_submit_calcfunction"
    nt.submit(wait=True)
    # print("results: ", results[])
    assert nt.nodes["sumdiff2"].node.outputs.sum == 9
