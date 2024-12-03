from IPython.display import IFrame


def test_workgraph_widget(wg_calcfunction):
    """Save the workgraph"""

    wg = wg_calcfunction
    wg.name = "test_workgraph_widget"
    wg.tasks["sumdiff2"].waiting_on.add(wg.tasks["sumdiff2"])
    value = wg.to_widget_value()
    assert len(value["nodes"]) == 2
    # the waiting_on is also transformed to links
    assert len(value["links"]) == 2
    # to_html
    data = wg.to_html()
    assert isinstance(data, IFrame)


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
