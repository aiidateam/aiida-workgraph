import pytest
import sys
import logging
import os

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main():
    """Run pytest with Playwright in non-headless mode."""
    pytest_args = [
        "tests/test_map.py::test_map_instruction",
        "-sv",  # Show output and verbose mode
    ]

    # Set the Playwright headless mode to 'no'
    os.environ["PYTEST_PLAYWRIGHT_HEADLESS"] = "no"

    logging.info(f"Running tests with arguments: {pytest_args}")

    # Run pytest with extracted arguments
    result = pytest.main(pytest_args)

    if result != 0:  # `pytest.main()` directly returns an integer exit code
        logging.error(f"Tests failed with exit code: {result}")
        sys.exit(result)
    else:
        logging.info("All tests passed successfully.")


if __name__ == "__main__":
    main()
