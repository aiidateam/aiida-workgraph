import aiida

aiida.load_profile()


def test_submit(wt_calcjob):
    """Submit simple calcjob."""
    wt = wt_calcjob
    wt.name = "test_submit_calcjob"
    wt.submit(wait=True)
    assert wt.nodes["add2"].outputs["sum"].value == 9
