# Sample test case for the root route
def test_root_route(client):
    response = client.get("/api")
    assert response.status_code == 200
    assert response.json() == {"message": "Welcome to AiiDA-WorkGraph."}


# Sample test case for the root route
def test_workgraph_route(client, wt_calcfunction):
    wt_calcfunction.submit(wait=True)
    response = client.get("/api/workgraph-data")
    assert response.status_code == 200
    assert len(response.json()) > 0
