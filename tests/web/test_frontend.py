import pytest

from playwright.sync_api import sync_playwright

import uvicorn

from multiprocessing import Process, Value

import contextlib
import threading
import time
import os

################################
# Utilities for frontend tests #
################################


class UvicornTestServer(uvicorn.Server):
    """
    Suggested way to start a server in a background by developers
    https://github.com/encode/uvicorn/discussions/1103#discussioncomment-941726
    """

    def install_signal_handlers(self):
        pass

    @contextlib.contextmanager
    def run_in_thread(self):
        thread = threading.Thread(target=self.run)
        thread.start()
        try:
            while not self.started:
                time.sleep(1e-3)
            yield
        finally:
            self.should_exit = True
            thread.join()


def run_uvicorn_web_server(
    web_server_started, stop_web_server, **uvicorn_configuration
):
    config = uvicorn.Config(**uvicorn_configuration)
    uvicorn_web_server = UvicornTestServer(config=config)
    with uvicorn_web_server.run_in_thread():
        with web_server_started.get_lock():
            web_server_started.value = 1
        while not stop_web_server.value:
            time.sleep(1e-3)


###############################
# Fixtuers for frontend tests #
###############################


@pytest.fixture(scope="module")
def uvicorn_configuration():
    return {
        "app": "aiida_workgraph.web.backend.app.api:app",
        "host": "0.0.0.0",
        "port": 8000,
        "log_level": "info",
        "workers": 2,
    }


@pytest.fixture(scope="module")
def web_server(set_backend_server_settings, uvicorn_configuration):
    from ctypes import c_int8

    web_server_started = Value(c_int8, 0)
    stop_web_server = Value(c_int8, 0)
    web_server_proc = Process(
        target=run_uvicorn_web_server,
        args=(web_server_started, stop_web_server),
        kwargs=uvicorn_configuration,
    )

    web_server_proc.start()

    while not web_server_started.value:
        time.sleep(1e-3)

    yield web_server_proc

    with stop_web_server.get_lock():
        stop_web_server.value = 1

    web_server_proc.join()
    web_server_proc.close()


# Define a fixture for the browser
@pytest.fixture(scope="module")
def browser():
    pytest_playwright_headless = os.environ.get("PYTEST_PLAYWRIGHT_HEADLESS", "yes")
    if pytest_playwright_headless == "yes":
        headless = True
    elif pytest_playwright_headless == "no":
        headless = False
    else:
        raise ValueError(
            f"Found environment variable PYTEST_PLAYWRIGHT_HEADLESS={pytest_playwright_headless}, "
            'please use "yes" or "no"'
        )

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        yield browser
        browser.close()


# Define a fixture for the page
@pytest.fixture(scope="module")
def page(browser):
    with browser.new_page() as page:
        yield page
        page.close()


##################
# Frontend tests #
##################


@pytest.mark.frontend
def test_homepage(web_server, page):
    page.goto("http://localhost:8000")

    assert page.title() == "AiiDA-WorkGraph App"

    # Check for the existence of a specific element on the page
    # Attempt to locate the element
    element = page.locator("a[href='/workgraph']")

    # Check if the element is found
    if not element.is_visible():
        pytest.fail("Element 'a[href='/wortre']' not found on the page")


@pytest.mark.frontend
def test_workgraph(web_server, page, aiida_profile, wg_calcfunction):
    wg_calcfunction.run()

    page.goto("http://localhost:8000")
    # Since the routing is done by react-router-dom we cannot access it with a call like this
    # page.goto("http://localhost:8000/workgraph" but have to navigate to it
    page.click('a[href="/workgraph"]')

    # Check for the existence of a specific element on the page

    # Verify the presence of the WorkGraphTable heading
    assert page.locator("h2").inner_text() == "WorkGraph"

    # Verify the presence of the search input
    assert page.locator(".search-input").is_visible()

    # Verify the presence of the table header columns
    # Verify the presence of the table header columns
    assert page.locator("th:has-text('PK')").is_visible()
    assert page.locator("th:has-text('Created')").is_visible()
    assert page.locator("th:has-text('Process Label')").is_visible()
    assert page.locator("th:has-text('State')").is_visible()
    assert page.locator("th:has-text('Actions')").is_visible()

    # Verify the presence of pagination controls
    assert page.locator(".pagination").is_visible()

    # Verify the presence of at least one row in the table
    page.wait_for_timeout(8000)
    assert page.locator("tr").count() >= 1  # Including header row


@pytest.mark.frontend
def test_workgraph_item(page, wg_calcfunction):
    wg = wg_calcfunction
    wg.run()
    page.goto("http://localhost:8000/workgraph/{}".format(wg.pk))
    page.wait_for_timeout(8000)

    page.get_by_text("sumdiff3").is_visible()

    # Simulate user interaction (e.g., clicking a button)
    # Replace the selector with the actual selector of the button you want to click
    # You should identify the button that triggers an action in your component
    page.get_by_role("button", name="Arrange").click()
    page.wait_for_timeout(8000)
    # Capture a screenshot
    screenshot = page.screenshot()

    # Save the screenshot to a file
    with open("screenshot.png", "wb") as f:
        f.write(screenshot)
