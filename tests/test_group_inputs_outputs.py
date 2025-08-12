from aiida_workgraph import WorkGraph, task


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


def test_load_from_db():
    """Test loading a WorkGraph from the database."""
    from aiida_workgraph.tasks.builtins import GraphInputs

    wg = WorkGraph("test_load_from_db")
    wg.inputs = {"x": 1, "y": 2}
    wg.save()
    wg2 = WorkGraph.load(wg.pk)
    wg2.restart()
    assert isinstance(wg2.tasks.graph_inputs, GraphInputs)
    wg2.inputs.z = 3
    wg2.save()
    wg3 = WorkGraph.load(wg2.pk)
    assert wg3.inputs.x == 1
    assert wg3.inputs.z == 3


def test_detect_graph_inputs(decorated_add):
    """Test that graph inputs are detected correctly."""

    # case 1: multiple graph inputs share the same value
    @task.graph
    def graph1(x, y):
        decorated_add(x=x, y=y)

    wg = graph1.build_graph(x=1, y=1)
    assert "graph_inputs.x -> add.x" in wg.links
    assert "graph_inputs.y -> add.y" in wg.links

    # case 2: input variable was renamed
    @task.graph
    def graph1(x, y):
        z = y
        decorated_add(x=x, y=z)

    assert "graph_inputs.x -> add.x" in wg.links
    assert "graph_inputs.y -> add.y" in wg.links
