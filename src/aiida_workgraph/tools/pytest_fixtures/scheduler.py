"""Fixtures to interact with the scheduler."""

from __future__ import annotations

import logging
import typing as t

import pytest

if t.TYPE_CHECKING:
    from aiida_workgraph.engine.scheduler.client import SchedulerClient

DEFAULT_TEST_SCHEDULER_NAME = "pytest-scheduler"


@pytest.fixture(scope="session")
def scheduler_client(
    aiida_profile, scheduler_name=DEFAULT_TEST_SCHEDULER_NAME
) -> "SchedulerClient":
    """Return a scheduler client for the configured test profile for the test session.

    The scheduler will be automatically stopped at the end of the test session.

    Usage::

        def test(scheduler_client):
            from aiida_workgraph.engine.scheduler.client import SchedulerClient
            assert isinstance(scheduler_client, SchedulerClient)

    """
    from aiida_workgraph.engine.scheduler import get_scheduler_client
    from aiida_workgraph.engine.scheduler.client import get_scheduler_node
    from aiida.engine.daemon.client import (
        DaemonNotRunningException,
        DaemonTimeoutException,
    )

    scheduler = get_scheduler_node(name=scheduler_name)

    if scheduler:
        raise ValueError(
            f"Scheduler {scheduler_name} is already existing. "
            f"Please remove it before running the test suite."
        )

    scheduler_client = get_scheduler_client(
        scheduler_name=scheduler_name, profile_name=aiida_profile.name
    )

    try:
        yield scheduler_client
    finally:
        try:
            scheduler_client.stop_daemon(wait=True)
        except DaemonNotRunningException:
            pass
        # Give an additional grace period by manually waiting for the scheduler to be stopped. In certain unit test
        # scenarios, the built in wait time in ``scheduler_client.stop_daemon`` is not sufficient and even though the
        # scheduler is stopped, ``scheduler_client.is_daemon_running`` will return false for a little bit longer.
        scheduler_client._await_condition(
            lambda: not scheduler_client.is_daemon_running,
            DaemonTimeoutException("The scheduler failed to stop."),
        )


@pytest.fixture
def started_scheduler_client(scheduler_client: "SchedulerClient"):
    """Ensure that the scheduler is running for the test profile and return the associated client.

    Usage::

        def test(started_scheduler_client):
            assert started_scheduler_client.is_daemon_running

    """
    from aiida_workgraph.engine.scheduler.client import (
        get_scheduler_node,
    )
    import time

    if not scheduler_client.is_daemon_running:
        scheduler_client.start_daemon()
        assert scheduler_client.is_daemon_running

    logger = logging.getLogger("tests.scheduler:started_scheduler_client")
    logger.debug(f"Daemon log file is located at: {scheduler_client.daemon_log_file}")
    # wait for scheduler to be started
    scheduler = get_scheduler_node(name=DEFAULT_TEST_SCHEDULER_NAME)
    timeout = 10
    tstart = time.time()
    while not scheduler:
        scheduler = get_scheduler_node(name=DEFAULT_TEST_SCHEDULER_NAME)
        time.sleep(0.5)
        if time.time() - tstart > timeout:
            raise TimeoutError(f"Scheduler failed to start after {timeout} seconds.")

    yield scheduler_client


@pytest.fixture
def stopped_scheduler_client(scheduler_client: "SchedulerClient"):
    """Ensure that the scheduler is not running for the test profile and return the associated client.

    Usage::

        def test(stopped_scheduler_client):
            assert not stopped_scheduler_client.is_daemon_running

    """
    from aiida_workgraph.engine.scheduler.client import DaemonTimeoutException

    if scheduler_client.is_daemon_running:
        scheduler_client.stop_daemon(wait=True)
        # Give an additional grace period by manually waiting for the scheduler to be stopped. In certain unit test
        # scenarios, the built in wait time in ``scheduler_client.stop_daemon`` is not sufficient and even though the
        # scheduler is stopped, ``scheduler_client.is_daemon_running`` will return false for a little bit longer.
        scheduler_client._await_condition(
            lambda: not scheduler_client.is_daemon_running,
            DaemonTimeoutException("The scheduler failed to stop."),
        )

    yield scheduler_client
