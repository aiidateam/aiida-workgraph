import aiida

aiida.load_profile()


def test_AtomsData():
    from aiida_workgraph.orm.atoms import AtomsData
    from ase.build import bulk

    atoms = bulk("Si")
    data = AtomsData(atoms)
    data.store()
    assert data.value == atoms
    assert data.base.attributes.get("formula") == "Si2"
    assert data.base.attributes.get("pbc") == [True, True, True]
