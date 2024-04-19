from aiida_workgraph import WorkGraph
import aiida

aiida.load_profile()


def test_numpy_array(decorated_normal_add):
    """Test data type with numpy array."""
    import numpy as np

    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    wg = WorkGraph()
    wg.nodes.new(decorated_normal_add, name="add1", x=x, y=y)
    wg.submit(wait=True)
    # wg.run()
    assert wg.state.upper() == "FINISHED"
    assert (wg.nodes["add1"].outputs["result"].value == np.array([5, 7, 9])).all()
