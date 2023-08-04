import aiida

aiida.load_profile()


def test_run(wt_calcfunction):
    """Run simple calcfunction."""
    wt = wt_calcfunction
    wt.name = "test_run_calcfunction"
    wt.run()
    print("state: ", wt.state)
    # print("results: ", results[])
    assert wt.nodes["sumdiff2"].node.outputs.sum == 9


def test_submit(wt_calcfunction):
    """Submit simple calcfunction."""
    wt = wt_calcfunction
    wt.name = "test_submit_calcfunction"
    wt.submit(wait=True)
    # print("results: ", results[])
    assert wt.nodes["sumdiff2"].node.outputs.sum == 9
