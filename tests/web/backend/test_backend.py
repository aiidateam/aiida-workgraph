import pytest


@pytest.mark.backend
def test_root_route(client):
    """Sample test case for the root route"""
    response = client.get("/api")
    assert response.status_code == 200
    assert response.json() == {"message": "Welcome to AiiDA-WorkGraph."}


@pytest.mark.backend
def test_workgraph_route(client, wg_calcfunction):
    """Sample test case for the root route"""
    wg_calcfunction.run()
    response = client.get("/api/workgraph-data")
    assert response.status_code == 200
    assert len(response.json()) > 0
