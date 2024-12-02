def test_workgraph_widget(wg_calcfunction):
    """Save the workgraph"""
    from IPython.display import IFrame

    wg = wg_calcfunction
    wg.name = "test_workgraph_widget"
    wg.tasks["sumdiff2"].waiting_on.add(wg.tasks["sumdiff2"])
    wg._widget.from_workgraph(wg)
    assert len(wg._widget.value["nodes"]) == 2
    # the waiting_on is also transformed to links
    assert len(wg._widget.value["links"]) == 2
    # to_html
    data = wg._widget.to_html()
    assert isinstance(data, IFrame)


def test_workgraph_task(wg_calcfunction):
    """Save the workgraph"""
    wg = wg_calcfunction
    wg.name = "test_workgraph_task"
    wg.tasks["sumdiff2"]._widget.from_node(wg.tasks["sumdiff2"])
    print(wg.tasks["sumdiff2"]._widget.value)
    assert len(wg.tasks["sumdiff2"]._widget.value["nodes"]) == 1
    assert len(
        wg.tasks["sumdiff2"]._widget.value["nodes"]["sumdiff2"]["inputs"]
    ) == len(wg.tasks["sumdiff2"].inputs)
    assert len(wg.tasks["sumdiff2"]._widget.value["links"]) == 0
