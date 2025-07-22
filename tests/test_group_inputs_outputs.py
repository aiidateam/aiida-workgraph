from aiida_workgraph import WorkGraph


def test_group_inputs_outputs(decorated_add):
    """Group inputs and outputs in a WorkGraph."""
    wg = WorkGraph("test_group_inputs_outputs")
    wg.inputs = {
        "add": {
            "x": 1,
            "y": 2,
        },
    }
    task1 = wg.add_task(decorated_add, x=wg.inputs.add.x, y=wg.inputs.add.y)
    task2 = wg.add_task(decorated_add, x=wg.inputs.add.x, y=task1.outputs.result)
    wg.outputs.results = {
        "sum1": task1.outputs.result,
        "sum2": task2.outputs.result,
    }
    wg.run()
    assert wg.outputs.results.sum1 == 4
    assert wg.outputs.results.sum2 == 6
    # the graph inputs will be serialized as AiiDA nodes
    assert (
        wg.process.inputs.workgraph_data.tasks.graph_inputs.inputs.sockets.add.sockets.x.property.value
        == 1
    )
