import aiida

aiida.load_profile()


def test_submit(wg_calcjob):
    """Submit simple calcjob."""
    wg = wg_calcjob
    wg.name = "test_submit_calcjob"
    wg.submit(wait=True)
    assert wg.nodes["add2"].outputs["sum"].value == 9
