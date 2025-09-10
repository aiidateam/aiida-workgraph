from aiida import orm


def test_get_or_create_code(fixture_localhost):
    from aiida_workgraph.utils import get_or_create_code
    from aiida.orm import Code

    # create a new code
    code1 = get_or_create_code(
        computer="localhost",
        code_label="test_code",
        code_path="/bin/bash",
        prepend_text='echo "Hello, World!"',
    )
    assert isinstance(code1, Code)
    # use already created code
    code2 = get_or_create_code(
        computer="localhost",
        code_label="test_code",
        code_path="/bin/bash",
        prepend_text='echo "Hello, World!"',
    )
    assert code1.uuid == code2.uuid


def test_get_parent_workgraphs():
    from aiida.common.links import LinkType
    from aiida_workgraph.utils import get_parent_workgraphs

    wn1 = orm.WorkflowNode()
    wn2 = orm.WorkflowNode()
    wn3 = orm.WorkflowNode()
    wn3.base.links.add_incoming(wn2, link_type=LinkType.CALL_WORK, link_label="link")
    wn2.base.links.add_incoming(wn1, link_type=LinkType.CALL_WORK, link_label="link")
    wn1.store()
    wn2.store()
    wn3.store()

    parent_workgraphs = get_parent_workgraphs(wn3.pk)
    assert len(parent_workgraphs) == 3


def test_generate_provenance_graph():

    from IPython.display import IFrame
    import os

    wn1 = orm.WorkflowNode()
    wn1.store()

    graph = generate_provenance_graph(wn1.pk)
    assert isinstance(graph, IFrame)
    # check file html/node_graph_{pk}.html is created
    assert os.path.isfile(f"html/node_graph_{wn1.pk}.html")
