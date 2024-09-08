import pytest
import re

from playwright.sync_api import expect


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
def test_workgraph(web_server, page, ran_wg_calcfunction):
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

    # Ensures that the first row has appeared
    page.get_by_role("cell", name="WorkGraph<test_debug_math>").hover()

    header_and_rows = page.get_by_role("row").all()
    # we wait for the cell to appear
    assert len(header_and_rows) == 2


@pytest.mark.frontend
def test_workgraph_item(web_server, page, ran_wg_calcfunction):
    page.goto("http://localhost:8000/workgraph/")
    page.get_by_role("link", name=str(ran_wg_calcfunction.pk), exact=True).click()
    # page.goto("http://localhost:8000/workgraph/{}".format(ran_wg_calcfunction.pk))
    # page.get_by_text()
    # page.get_by_role("button", name="Arrange").click()
    # ran_wg_calcfunction.pk))

    page.get_by_text("sumdiff3").is_visible()

    # Simulate user interaction (e.g., clicking a button)
    # Replace the selector with the actual selector of the button you want to click
    # You should identify the button that triggers an action in your component
    page.get_by_role("button", name="Arrange").click()

    gui_node = page.get_by_text("sumdiff2")

    # Verify that clicking on the "Real-time state" changes the color to green
    gui_node_color = gui_node.evaluate(
        "element => window.getComputedStyle(element).backgroundColor"
    )
    assert gui_node_color == "rgba(0, 0, 0, 0)"
    page.locator(".realtime-switch").click()

    # this waits until a green background appears
    page.wait_for_function(
        "selector => !!document.querySelector(selector)",
        arg="div.title[style='background: green;']",
    )
    gui_node_color = gui_node.evaluate(
        "element => window.getComputedStyle(element).backgroundColor"
    )
    assert gui_node_color == "rgb(0, 128, 0)"

    # Check the node-detail-view switch
    page.locator(".detail-switch").click()
    # pause here for debugging
    # page.pause()
    # page.wait_for_selector('[data-testid="input-x"] input', timeout=5000)  # Wait for up to 5 seconds
    # input_x_control = page.get_by_test_id("input-x").first.locator("input")
    # assert input_x_control.input_value() == "2"

    # Verify that clicking on the gui node will pop up a sidebar
    gui_node.click()
    node_details_sidebar = page.get_by_text("CloseNode")
    node_details_sidebar.wait_for(state="visible")
    assert node_details_sidebar.is_visible()
    node_details_sidebar.get_by_role("button", name="Close").click()
    node_details_sidebar.wait_for(state="hidden")
    assert node_details_sidebar.is_hidden()

    # verify Summary works
    page.get_by_role("button", name="Summary").click()
    assert page.get_by_text("typeWorkGraph<test_debug_math>").is_visible()

    # Verify that Log  works
    page.get_by_role("button", name="Log").click()
    log_line = (
        page.locator(".log-content")
        .locator("div")
        .filter(has_text=re.compile(r".*Task: sumdiff2 finished.*"))
    )
    log_line.wait_for(state="visible")
    assert log_line.is_visible()

    # Verify that Time  works
    page.get_by_role("button", name="Time").click()
    row = page.locator(".rct-sidebar-row ").get_by_text("sumdiff2")
    row.wait_for(state="visible")
    assert row.is_visible()


@pytest.mark.frontend
def test_datanode_item(web_server, page, ran_wg_calcfunction):
    page.goto("http://localhost:8000/datanode/")
    data_node_pk = ran_wg_calcfunction.nodes["sumdiff1"].inputs["x"].value.pk
    page.get_by_role("link", name=str(data_node_pk), exact=True).click()

    # check if three rows (header plus 2) are present
    expect(page.locator(":nth-match(tr, 3)")).to_be_visible()
    rows = page.get_by_role("row").all()
    assert "value" in rows[1].text_content()
    assert "node_type" in rows[2].text_content()


@pytest.mark.frontend
def test_settings(web_server, page, ran_wg_calcfunction):
    page.goto("http://localhost:8000/settings/")
    # Verify that only one row is visible
    expect(page.locator(":nth-match(tr, 1)")).to_be_visible()
    expect(page.locator(":nth-match(tr, 2)")).to_be_hidden()

    # Verify that after starting the daemon one additional row appeared
    page.get_by_role("button", name="Start Daemon").click()
    expect(page.locator(":nth-match(tr, 2)")).to_be_visible()
    expect(page.locator(":nth-match(tr, 3)")).to_be_hidden()

    # Verify that after adding workers one additional row appeared
    page.get_by_role("button", name="Increase Workers").click()
    expect(page.locator(":nth-match(tr, 3)")).to_be_visible()
    expect(page.locator(":nth-match(tr, 4)")).to_be_hidden()

    # Verify that after decreasing workers one row disappears
    page.get_by_role("button", name="Decrease Workers").click()
    expect(page.locator(":nth-match(tr, 2)")).to_be_visible()
    expect(page.locator(":nth-match(tr, 3)")).to_be_hidden()

    # Verify that stopping the daemon only the header row exists
    page.get_by_role("button", name="Stop Daemon").click()
    expect(page.locator(":nth-match(tr, 1)")).to_be_visible()
    expect(page.locator(":nth-match(tr, 2)")).to_be_hidden()


############################
# Tests mutating the state #
############################
# The tests below change the state of the aiida database and cannot necessary be executed in arbitrary order.


@pytest.mark.frontend
def test_workgraph_delete(web_server, page, ran_wg_calcfunction):
    """Tests that the last workgraph node in the table can be deleted successfully."""
    page.goto("http://localhost:8000")
    # Since the routing is done by react-router-dom we cannot access it with a call like this
    # page.goto("http://localhost:8000/workgraph" but have to navigate to it
    page.click('a[href="/workgraph"]')

    # Ensures that the last data row has appeared, the first row is header
    last_row = page.locator(":nth-match(tr, 2)")
    expect(last_row).to_be_visible()
    # verify that this is the last row
    expect(page.locator(":nth-match(tr, 3)")).to_be_hidden()

    delete_button = last_row.locator(".delete-button")

    # Verify that cancel or closing the prompt does not delete the node
    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_label("Close").click()
    expect(last_row).to_be_visible()

    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_role("button", name="Cancel").click()
    expect(last_row).to_be_visible()

    # Verify that confirming the prompt does delete the node
    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_role("button", name="Confirm").click()
    expect(last_row).to_be_hidden()


@pytest.mark.frontend
def test_datanode_delete(web_server, page, ran_wg_calcfunction):
    """Tests that the last data node in the table can be deleted successfully."""
    page.goto("http://localhost:8000")
    # Since the routing is done by react-router-dom we cannot access it with a call like this
    # page.goto("http://localhost:8000/workgraph" but have to navigate to it
    page.click('a[href="/datanode"]')

    # Ensures that the last data row has appeared
    last_row = page.locator(":nth-match(tr, 7)")
    expect(last_row).to_be_visible()
    # verify that this is the last row
    expect(page.locator(":nth-match(tr, 8)")).to_be_hidden()

    delete_button = last_row.locator(".delete-button")

    # Verify that cancel or closing the prompt does not delete the node
    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_label("Close").click()
    expect(last_row).to_be_visible()

    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_role("button", name="Cancel").click()
    expect(last_row).to_be_visible()

    # Verify that confirming the prompt does delete the node
    delete_button.click()
    expect(page.get_by_text("Confirm deletion")).to_be_visible()
    page.get_by_role("button", name="Confirm").click()
    expect(last_row).to_be_hidden()
