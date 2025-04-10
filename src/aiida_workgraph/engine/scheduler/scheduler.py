from __future__ import annotations

import asyncio
import functools
import logging
import threading
import uuid
from typing import Any, Dict, Optional, Union

import kiwipy
from plumpy.communications import wrap_communicator
from plumpy.events import set_event_loop_policy
from plumpy.process_comms import (
    RemoteProcessThreadController,
    TASK_ARGS,
    TASK_KEY,
    CONTINUE_TASK,
)
from aiida.common import exceptions
from aiida.orm import load_node
from aiida.engine.processes import ProcessState
from aiida.manage import get_manager
from aiida.brokers.rabbitmq.utils import (
    get_message_exchange_name,
    get_task_exchange_name,
)
from aiida_workgraph.orm.scheduler import SchedulerNode

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

LOGGER = logging.getLogger(__name__)

__all__ = ("Scheduler",)

SCHEDULER_PRIORITY_KEY = "_scheduler_priority"


class Scheduler:
    """
    A Scheduler class that:
      - Receives 'continue' tasks via RabbitMQ,
      - Runs them up to a configured concurrency limit,
      - Persists its queue/active-state in a SchedulerNode,
      - Assigns priorities to processes and always launches the
        highest-priority (lowest integer) first.
    """

    _communicator: Optional[kiwipy.Communicator] = None
    _controller: Optional[RemoteProcessThreadController] = None

    def __init__(
        self,
        name: str,
        maxium_calcjob: Optional[int] = None,
        poll_interval: Union[int, float] = 0,
        reset: bool = False,
    ):
        """
        :param name: A unique name for this scheduler.
        :param maxium_calcjob: Maximum number of processes to run concurrently
        :param poll_interval: Interval in seconds for the fallback polling mechanism
                              (if an AiiDA broadcast is missed).
        """
        self.name = name
        self._node: Optional[SchedulerNode] = None
        self._poll_interval = poll_interval

        # Set up AiiDA profile and event loop policy
        manager = get_manager()
        self._profile = manager.get_profile()
        self._prefix = f"aiida-{self._profile.uuid}"
        set_event_loop_policy()

        # We create our dedicated event loop
        self._loop = asyncio.new_event_loop()

        # Set up the RabbitMQ communicator and attach the task receiver callback
        self.communicator.add_task_subscriber(self.task_receiver)
        self.process_controller = get_manager().get_process_controller()

        if reset:
            self.reset()

        # If concurrency limit was provided, store it
        if maxium_calcjob is not None:
            self.node.maxium_calcjob = maxium_calcjob

    @property
    def node(self) -> SchedulerNode:
        """
        Return (and lazily create) the SchedulerNode that stores all data for this scheduler:
          - waiting_process: list of PKs waiting to run
          - running_process: list of PKs currently running
          - maxium_calcjob: concurrency limit
          - next_priority: the integer to assign to the next new process (defaults to 0)
        """
        if self._node is None:
            from aiida.orm import QueryBuilder

            qb = QueryBuilder()
            qb.append(SchedulerNode, filters={"attributes.name": self.name})
            if qb.count() == 0:
                # Create a new SchedulerNode
                self._node = SchedulerNode()
                self._node.name = self.name
                self._node.next_priority = 0
                self._node.store()
                LOGGER.info(
                    f"Created new SchedulerNode for {self.name}, pk={self._node.pk}"
                )
            else:
                # Load existing SchedulerNode
                self._node = qb.first()[0]
                LOGGER.info(
                    "Loaded existing SchedulerNode<%d> for '%s'",
                    self._node.pk,
                    self.name,
                )
                # Check for any running processes that might have finished
                self._check_for_finished_running()

        return self._node

    def reset(self) -> None:
        """
        Reset the scheduler by clearing the waiting and running process lists.
        """
        self.node.waiting_process = []
        self.node.running_process = []
        self.node.next_priority = 0
        LOGGER.info("Scheduler '%s' reset.", self.name)

    def _get_next_priority(self) -> int:
        """
        Return and then increment the 'next_priority' counter
        in the SchedulerNode’s extras or attributes.
        """
        current = self._node.next_priority
        self._node.next_priority = current - 1
        return current

    def get_url(self) -> str:
        """Return the RabbitMQ URL for this AiiDA profile."""
        from aiida.brokers.rabbitmq.utils import get_rmq_url

        kwargs = {
            key[7:]: val
            for key, val in self._profile.process_control_config.items()
            if key.startswith("broker_")
        }
        additional_kwargs = kwargs.pop("parameters", {})
        return get_rmq_url(**kwargs, **additional_kwargs)

    @property
    def communicator(self) -> kiwipy.Communicator:
        """
        Access (and lazily create) a kiwipy Communicator connected to RabbitMQ.
        We wrap it in `wrap_communicator` so it is driven by our event loop.
        """
        from kiwipy.rmq import RmqThreadCommunicator
        from aiida.orm.utils import serialize

        if self._communicator is None:
            communicator = RmqThreadCommunicator.connect(
                connection_params={"url": self.get_url()},
                message_exchange=get_message_exchange_name(self._prefix),
                encoder=functools.partial(serialize.serialize, encoding="utf-8"),
                decoder=serialize.deserialize_unsafe,
                task_exchange=get_task_exchange_name(self._prefix),
                task_queue=self.name,
                task_prefetch_count=10000,
            )
            self._communicator = wrap_communicator(communicator, self._loop)
            self._controller = RemoteProcessThreadController(communicator)
        return self._communicator

    @property
    def controller(self) -> RemoteProcessThreadController:
        """
        Access the remote process thread controller, used by AiiDA to continue processes.
        """
        if self._controller is None:
            raise exceptions.InternalError(
                "Controller not set (communicator not yet initialized)"
            )
        return self._controller

    async def task_receiver(
        self, communicator: kiwipy.Communicator, task: Dict[str, Any]
    ) -> None:
        """
        Subscriber callback for tasks posted to our queue.
        We only handle CONTINUE_TASK tasks.
        """
        task_type = task[TASK_KEY]
        if task_type == CONTINUE_TASK:
            pid = task.get(TASK_ARGS, {}).get("pid")
            if pid is None:
                return

            LOGGER.info("Received CONTINUE_TASK for pk=%d", pid)
            priority = self._compute_priority_for_new_process(pid)
            child_node = load_node(pid)
            child_node.base.extras.set(SCHEDULER_PRIORITY_KEY, priority)
            self.node.append_waiting_process(pid)
            # Trigger consumption (start more processes if possible)
            self.consume_process_queue()

    def _compute_priority_for_new_process(self, pid: int) -> int:
        """"""
        from aiida.common.links import LinkType

        # find the parent process
        node = load_node(pid)
        links = node.base.links.get_incoming()
        # find the parent node
        parent_nodes = [
            link.node
            for link in links
            if link.link_type in [LinkType.CALL_CALC, LinkType.CALL_WORK]
        ]
        if parent_nodes:
            parent_node = parent_nodes[0]
            return parent_node.base.extras.get(SCHEDULER_PRIORITY_KEY, 0)
        return self._get_next_priority()

    def consume_process_queue(self) -> None:
        """
        Attempt to start processes from the waiting queue,
        up to the concurrency limit. We pick the next process
        with the *highest* priority.
        """
        waiting = len(self.node.waiting_process)
        running = len(self.node.running_calcjob)
        max_ = self.node.maxium_calcjob
        LOGGER.info(
            "Summary: waiting_process=%d, running_calcjob=%d, max_calcjob=%d",
            waiting,
            running,
            max_,
        )

        # While we still have capacity, pop from waiting and continue
        while len(self.node.running_calcjob) < self.node.maxium_calcjob:
            # pick the next waiting PK with the highest priority
            next_pk = self._pop_highest_priority_waiting_process()
            if next_pk is None:
                LOGGER.info("No more processes in waiting queue.")
                return
            self.continue_process(next_pk)

        LOGGER.info(
            "Maximum concurrency (%d) reached, waiting for a calcjob to finish...",
            self.node.maxium_calcjob,
        )

    def _pop_highest_priority_waiting_process(self) -> Optional[int]:
        """Return the waiting process PK with the largest 'priority' extras, or None if empty."""
        waiting_list = self.node.waiting_process
        if not waiting_list:
            return None

        # Sort waiting PKs by their priority
        # The process with the largest integer is considered highest priority
        priority_list = [self._get_priority_for_pid(pk) for pk in waiting_list]
        # for i, pk in enumerate(waiting_list):
        # LOGGER.info(f"Waiting process pk={pk} has priority {priority_list[i]}")
        best_pk = waiting_list[priority_list.index(max(priority_list))]
        LOGGER.info(f"Best waiting process pk={best_pk}")
        self.node.remove_waiting_process(best_pk)
        return best_pk

    def _get_priority_for_pid(self, pk: int) -> int:
        """Return the integer priority from the node extras, or a large default if missing."""
        child_node = load_node(pk)
        return child_node.base.extras.get(SCHEDULER_PRIORITY_KEY, 9999999)

    def _check_for_finished_running(self) -> None:
        """
        If the scheduler was stopped while processes were running,
        some might have already finished. We'll remove them from 'running'
        if they're terminated. We also re-attach callbacks to any that are still running.
        """
        to_remove = []
        for pk in list(self.node.running_process):
            try:
                node = load_node(pk)
                if node.is_terminated:
                    LOGGER.info(
                        "Found pk=%d in 'running' but it's already terminated. Removing it.",
                        pk,
                    )
                    to_remove.append(pk)
                else:
                    # Re-attach callback in case we never did or the scheduler was restarted
                    self.call_on_process_finish(pk)
            except exceptions.NotExistent:
                LOGGER.warning(
                    "Found pk=%d in 'running' but it doesn't exist anymore. Removing it.",
                    pk,
                )
                to_remove.append(pk)

        for pk in to_remove:
            self.node.remove_running_process(pk)

    def continue_process(self, pk: int) -> None:
        """
        Actually continue an AiiDA process in the daemon. Attach a callback
        so that when it finishes, we free up a slot and can launch a new one.
        """
        LOGGER.info("Continuing process pk=%d...", pk)
        self.call_on_process_finish(pk)
        self.node.append_running_process(pk)
        self.process_controller.continue_process(pk)

    def call_on_process_finish(self, pk: int) -> None:
        """
        Attach both a broadcast-based subscriber and a fallback polling
        so that when process <pk> terminates, we call `on_process_finished_callback`.
        """
        node = load_node(pk=pk)
        subscriber_identifier = str(uuid.uuid4())
        done_event = threading.Event()

        # Our final callback is scheduled in the loop so it runs in the correct thread
        callback = functools.partial(
            self._loop.call_soon, self.on_process_finished_callback, pk
        )

        def inline_callback(event, *args, **kwargs):
            if event.is_set():
                return
            try:
                callback()
            finally:
                event.set()
                if self.communicator:
                    self.communicator.remove_broadcast_subscriber(subscriber_identifier)

        # Subscribe to broadcast events from that PK
        broadcast_filter = kiwipy.BroadcastFilter(
            functools.partial(inline_callback, done_event), sender=pk
        )
        for state in [
            ProcessState.FINISHED,
            ProcessState.KILLED,
            ProcessState.EXCEPTED,
        ]:
            broadcast_filter.add_subject_filter(f"state_changed.*.{state.value}")

        # Register the broadcast subscriber
        if self.communicator:
            LOGGER.debug("Adding broadcast subscriber for pk=%d", pk)
            self.communicator.add_broadcast_subscriber(
                broadcast_filter, subscriber_identifier
            )

        # Start polling in parallel, as a fallback
        self._poll_process(node, functools.partial(inline_callback, done_event))

    def on_process_finished_callback(self, pk: int) -> None:
        """
        When a process finishes, remove it from the running list, and try
        to consume another from the waiting queue.
        """
        LOGGER.info("[📢] Process pk=%d finished", pk)
        if pk in self.node.running_process:
            self.node.remove_running_process(pk)
        # Attempt to start more processes if any remain
        self.consume_process_queue()

    def _poll_process(self, node, callback):
        """
        Fallback: check if the process is terminated; if so call `callback`,
        otherwise schedule another check in `self._poll_interval` seconds.
        """
        if node.is_terminated:
            LOGGER.debug(
                "%s<%d> is already terminated (confirmed by polling).",
                node.__class__.__name__,
                node.pk,
            )
            self._loop.call_soon(callback)
        else:
            self._loop.call_later(
                self._poll_interval, self._poll_process, node, callback
            )

    def run(self) -> None:
        """
        Single entry point for the user: blocks the current thread
        and runs the event loop forever, listening for tasks and
        automatically continuing any leftover processes from before.
        """
        LOGGER.info(
            "Starting Scheduler '%s' main loop. Press Ctrl+C to stop.", self.name
        )
        # Once loop starts, re-check leftover waiting/running
        self._loop.call_soon(self.consume_process_queue)

        try:
            self._loop.run_forever()
        except KeyboardInterrupt:
            LOGGER.info("Scheduler '%s' interrupted by user (Ctrl+C).", self.name)
        finally:
            self._loop.run_until_complete(self._loop.shutdown_asyncgens())
            self._loop.close()
            LOGGER.info("Scheduler '%s' event loop closed.", self.name)
