from typing import Dict
from aiida_workgraph import task
from aiida.orm import StructureData, UpfData


@task(
    inputs=[["String", "pseudo_family"], [StructureData, "structure"]],
    outputs=[[UpfData, "Pseudo"]],
)
def get_pseudo_from_structure(
    pseudo_family: str, structure: StructureData
) -> Dict[str, UpfData]:
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
