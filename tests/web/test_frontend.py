import pytest


def test_homepage(page):
    page.goto("http://localhost:8000")

    assert page.title() == "AiiDA-WorkGraph App"

    # Check for the existence of a specific element on the page
    # Attempt to locate the element
    element = page.locator("a[href='/workgraph']")

    # Check if the element is found
    if not element.is_visible():
        pytest.fail("Element 'a[href='/wortre']' not found on the page")


def test_workgraph(page, wg_calcfunction):
    wg_calcfunction.submit(wait=True)
    page.goto("http://localhost:8000/workgraph")

    # Check for the existence of a specific element on the page

    # Verify the presence of the WorkGraphTable heading
    assert page.locator("h2").inner_text() == "WorkGraph"

    # Verify the presence of the search input
    assert page.locator(".search-input").is_visible()

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
    assert page.locator("tr").count() >= 2  # Including header row


def test_workgraph_item(page, wg_calcfunction, assert_snapshot):

    wg = wg_calcfunction
    wg.submit(wait=True)
    page.goto("http://localhost:8000/workgraph/{}".format(wg.pk))
    page.wait_for_timeout(8000)

    page.get_by_text("sumdiff3").is_visible()

    # Simulate user interaction (e.g., clicking a button)
    page.get_by_role("button", name="Arrange").click()
    page.wait_for_timeout(8000)

    # compare the screenshot
    assert_snapshot(page.screenshot(), "screenshot.png")
