from aiida_workgraph import node
from aiida.orm import StructureData, UpfData


@node(
    inputs=[["String", "pseudo_family"], [StructureData, "structure"]],
    outputs=[[UpfData, "Pseudo"]],
)
def get_pseudo_from_structure(pseudo_family, structure):
    """for input_namespace"""
    from aiida.orm import Group, QueryBuilder

    pseudo_group = (
        QueryBuilder().append(Group, filters={"label": pseudo_family}).one()[0]
    )
    elements = [kind.name for kind in structure.kinds]
    pseudos = {}
    for ele in elements:
        for n in pseudo_group.nodes:
            if ele == n.element:
                pseudos[ele] = n
    return {"Pseudo": pseudos}


if __name__ == "__main__":
    pass
