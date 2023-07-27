import time


def wait_nt(nt, timeout=50):
    """Wait for node tree to finish"""
    start = time.time()
    nt.update()
    while nt.state not in ("PAUSED", "FINISHED", "FAILED", "CANCELLED"):
        time.sleep(0.5)
        if time.time() - start > timeout:
            break
        nt.update()
