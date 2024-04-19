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


def test_workgraph(page, wt_calcfunction):
    wt_calcfunction.submit(wait=True)
    page.goto("http://localhost:8000/workgraph")

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
    assert page.locator("tr").count() >= 2  # Including header row


def test_workgraph_item(page, wt_calcfunction):

    wt = wt_calcfunction
    wt.submit(wait=True)
    page.goto("http://localhost:8000/workgraph/{}".format(wt.pk))
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
