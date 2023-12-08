from aiida import load_profile

load_profile()


def get_pseudo_from_structure(pseudo_family, structure):
    """for input_namespace"""
    from aiida.orm import Group, QueryBuilder

    pseudo_group = (
        QueryBuilder().append(Group, filters={"label": pseudo_family}).one()[0]
    )
    elements = [kind.name for kind in structure.kinds]
    pseudos = {}
    print("elements: ", elements)
    for ele in elements:
        for node in pseudo_group.nodes:
            if ele == node.element:
                pseudos[ele] = node
    print("pseudos: ", pseudos)
    return pseudos


if __name__ == "__main__":
    pass
