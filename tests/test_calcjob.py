import time
import numpy as np
import aiida

aiida.load_profile()


def test_submit(nt_calcjob):
    """Submit simple calcjob."""
    nt = nt_calcjob
    nt.name = "test_submit_calcjob"
    nt.submit(wait=True)
    # print("results: ", results[])
    assert nt.nodes["add2"].node.outputs.sum == 9
