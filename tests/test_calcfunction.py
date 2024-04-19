import aiida

aiida.load_profile()


def test_run(wg_calcfunction):
    """Run simple calcfunction."""
    wg = wg_calcfunction
    wg.name = "test_run_calcfunction"
    wg.run()
    print("state: ", wg.state)
    # print("results: ", results[])
    assert wg.nodes["sumdiff2"].node.outputs.sum == 9
    assert wg.nodes["sumdiff2"].outputs["sum"].value == 9


def test_submit(wg_calcfunction):
    """Submit simple calcfunction."""
    wg = wg_calcfunction
    wg.name = "test_submit_calcfunction"
    wg.submit(wait=True)
    # print("results: ", results[])
    assert wg.nodes["sumdiff2"].outputs["sum"].value == 9
