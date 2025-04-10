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


class Scheduler:
    """
    A Scheduler class that receives 'continue' tasks via RabbitMQ,
    runs them up to a configured concurrency limit, and persists
    its queue/active-state in a SchedulerNode in the AiiDA database.
    """

    _communicator: Optional[kiwipy.Communicator] = None
    _controller: Optional[RemoteProcessThreadController] = None

    def __init__(
        self,
        name: str,
        maxium_running_processes: Optional[int] = None,
        poll_interval: Union[int, float] = 0,
    ):
        """
        :param name: A unique name for this scheduler.
        :param maxium_running_processes: Maximum number of processes to run concurrently
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

        # We create our dedicated event loop and keep it in a private variable
        self._loop = asyncio.new_event_loop()

        # Set up the RabbitMQ communicator and attach the task receiver callback
        self.communicator.add_task_subscriber(self.task_receiver)
        self.process_controller = get_manager().get_process_controller()

        # If the user provided a concurrency limit, update the node
        if maxium_running_processes is not None:
            self.node.maxium_running_processes = maxium_running_processes

    @property
    def node(self) -> SchedulerNode:
        """
        Return the SchedulerNode that stores all persistent data for this scheduler:
          - waiting_process (list of PKs waiting to run)
          - running_processes (list of PKs currently running)
          - maxium_running_processes (int concurrency limit)
        """
        if self._node is None:
            from aiida.orm import QueryBuilder

            qb = QueryBuilder()
            qb.append(SchedulerNode, filters={"attributes.name": self.name})
            if qb.count() == 0:
                # create new SchedulerNode
                self._node = SchedulerNode()
                self._node.name = self.name
                self._node.store()
                LOGGER.info(
                    f"Created new SchedulerNode for {self.name}, pk={self._node.pk}"
                )
            else:
                # load existing SchedulerNode
                self._node = qb.first()[0]
                LOGGER.info(
                    "Loaded existing SchedulerNode<%d> for '%s'",
                    self._node.pk,
                    self.name,
                )
                # Ensure any "running" processes are still actually running
                self._check_for_finished_running()
        return self._node

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
        We wrap it in `wrap_communicator` so it is driven by our internal event loop.
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
            if pid is not None:
                LOGGER.info("Received CONTINUE_TASK for pk=%d", pid)
                self.node.append_waiting_process(pid)
                # Trigger consumption (start more processes if possible)
                self.consume_process_queue()

    def consume_process_queue(self) -> None:
        """
        Attempt to start processes from the waiting queue,
        up to the concurrency limit in maxium_running_processes.
        """
        waiting = len(self.node.waiting_process)
        running = len(self.node.running_processes)
        max_ = self.node.maxium_running_processes
        LOGGER.info(
            "consume_process_queue: waiting=%d, running=%d, max=%d",
            waiting,
            running,
            max_,
        )

        # While we still have capacity, pop from waiting and continue
        while len(self.node.running_processes) < self.node.maxium_running_processes:
            if self.node.waiting_process:
                pk = self.node.pop_waiting_process()
                self.continue_process(pk)
            else:
                LOGGER.info("No more processes in waiting queue.")
                return

        LOGGER.info(
            "Maximum concurrency (%d) reached, waiting for a process to finish...",
            self.node.maxium_running_processes,
        )

    def _check_for_finished_running(self) -> None:
        """
        If the scheduler was stopped while processes were running,
        some might have already finished by now. Let's check each quickly.
        If it's finished, we remove it from 'running_processes'.
        """
        to_remove = []
        for pk in list(self.node.running_processes):
            node = load_node(pk)
            if node.is_terminated:
                LOGGER.info(
                    "Found pk=%d in 'running' but it is actually terminated. Removing it.",
                    pk,
                )
                to_remove.append(pk)
            else:
                # add the callback to the process
                self.call_on_process_finish(pk)
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
        so that when the process <pk> terminates, we call `on_process_finished_callback`.
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

        # Subscribe to broadcast events from that PK,
        # for any terminal state (FINISHED, KILLED, EXCEPTED)
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
        if pk in self.node.running_processes:
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

        # As soon as the loop starts, let's see if we have leftover waiting/running
        # processes and try to continue them. We schedule this with call_soon
        # so it runs right after the loop is in run_forever().
        self._loop.call_soon(self.consume_process_queue)

        try:
            self._loop.run_forever()
        except KeyboardInterrupt:
            LOGGER.info("Scheduler '%s' interrupted by user (Ctrl+C).", self.name)
        finally:
            # Cleanup
            self._loop.run_until_complete(self._loop.shutdown_asyncgens())
            self._loop.close()
            LOGGER.info("Scheduler '%s' event loop closed.", self.name)
