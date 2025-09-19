from aiida_workgraph import WorkGraph


def test_zone_task(decorated_add):
    """Test the zone task."""

    wg = WorkGraph('test_zone')
    add1 = wg.add_task(decorated_add, name='add1', x=1, y=1)
    zone1 = wg.add_task('workgraph.zone', name='zone1')
    zone1.add_task(decorated_add, name='add2', x=1, y=1)
    zone1.add_task(decorated_add, name='add3', x=1, y=add1.outputs.result)
    wg.add_task(decorated_add, name='add4', x=1, y=wg.tasks.add2.outputs.result)
    wg.add_task(decorated_add, name='add5', x=1, y=wg.tasks.add3.outputs.result)
    connectivity = wg.build_connectivity()
    assert connectivity['zone']['add4']['input_tasks'] == ['zone1']
    assert connectivity['zone']['add5']['input_tasks'] == ['zone1']
