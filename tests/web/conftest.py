import pytest
from fastapi.testclient import TestClient
from aiida_workgraph.web.backend.app.api import app
from playwright.sync_api import sync_playwright


# Define a fixture for the FastAPI app client
@pytest.fixture(scope="module")
def client():
    return TestClient(app)


# Define a fixture for the browser
@pytest.fixture(scope="module")
def browser():
    with sync_playwright() as p:
        browser = p.chromium.launch()
        yield browser
        browser.close()


# Define a fixture for the page
@pytest.fixture(scope="module")
def page(browser):
    with browser.new_page() as page:
        yield page
        page.close()
