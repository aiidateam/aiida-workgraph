import pytest
from fastapi.testclient import TestClient

##############################
# Fixtures for backend tests #
##############################


@pytest.fixture(scope="module")
def client(set_backend_server_settings):
    from aiida_workgraph.web.backend.app.api import app

    return TestClient(app)


#################
# Backend tests #
#################

# Sample test case for the root route
@pytest.mark.backend
def test_root_route(client):
    response = client.get("/api")
    assert response.status_code == 200
    assert response.json() == {"message": "Welcome to AiiDA-WorkGraph."}


# Sample test case for the root route
@pytest.mark.backend
def test_workgraph_route(client, wg_calcfunction):
    wg_calcfunction.run()
    response = client.get("/api/workgraph-data")
    assert response.status_code == 200
    assert len(response.json()) > 0
