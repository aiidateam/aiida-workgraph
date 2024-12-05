from IPython.display import IFrame


def test_workgraph_widget(wg_calcfunction, decorated_add):
    """Save the workgraph"""

    wg = wg_calcfunction
    wg.name = "test_workgraph_widget"
    wg.add_task(decorated_add, "add1", x=1, y=3)
    wg.tasks["sumdiff2"].waiting_on.add(wg.tasks["sumdiff2"])
    value = wg.to_widget_value()
    assert len(value["nodes"]) == 3
    # the waiting_on is also transformed to links
    assert len(value["links"]) == 2
    # check required sockets
    # there are more than 2 inputs, but only 2 are required
    assert len(wg.tasks["add1"].inputs) > 2
    assert len(value["nodes"]["add1"]["inputs"]) == 2
    # to_html
    data = wg.to_html()
    assert isinstance(data, IFrame)
    # check _repr_mimebundle_ is working
    data = wg._repr_mimebundle_()


def test_workgraph_task(wg_calcfunction):
    """Save the workgraph"""
    wg = wg_calcfunction
    value = wg.tasks["sumdiff2"].to_widget_value()
    assert len(value["nodes"]) == 1
    assert len(value["nodes"]["sumdiff2"]["inputs"]) == len(wg.tasks["sumdiff2"].inputs)
    assert len(value["links"]) == 0
    # to html
    data = wg.tasks["sumdiff2"].to_html()
    assert isinstance(data, IFrame)
    # check _repr_mimebundle_ is working
    data = wg.tasks["sumdiff2"]._repr_mimebundle_()
