# Sample test case for the root route
def test_root_route(client):
    response = client.get("/api")
    assert response.status_code == 200
    assert response.json() == {"message": "Welcome to AiiDA-WorkTree."}


# Add more test cases for your other routes and features
# For example, testing authentication, database interactions, etc.
