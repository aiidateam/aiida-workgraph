from __future__ import annotations

import asyncio
import functools
import logging
import threading
import uuid
from typing import Any, Dict, Optional, Union

import kiwipy
from plumpy.communications import wrap_communicator
from plumpy.events import set_event_loop_policy, reset_event_loop_policy
from plumpy.process_comms import (
    RemoteProcessThreadController,
    TASK_ARGS,
    TASK_KEY,
    CONTINUE_TASK,
)
from aiida.common import exceptions
from aiida.orm import load_node, QueryBuilder, Node, CalcJobNode
from aiida.engine.processes import ProcessState
from aiida.manage import get_manager
from aiida.brokers.rabbitmq.utils import (
    get_message_exchange_name,
    get_task_exchange_name,
)
from aiida_workgraph.orm.scheduler import SchedulerNode
from aiida_workgraph.utils import (
    query_terminated_processes,
    query_existing_processes,
)

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

LOGGER = logging.getLogger(__name__)

__all__ = ("Scheduler",)

SCHEDULER_PRIORITY_KEY = "_scheduler_priority"
SCHEDULER_KEY = "_scheduler"
INTENT_KEY = "intent"
MESSAGE_TEXT_KEY = "message"


class Intent:
    """Intent constants for a process message"""

    STOP: str = "stop"
    CONTINUE: str = "continue"
    PLAY: str = "play"
    REPAIR: str = "repair"
    STATUS: str = "status"
    SET: str = "set"


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
        max_calcjobs: Optional[int] = None,
        max_workflows: Optional[int] = None,
        max_processes: Optional[int] = None,
        poll_interval: Union[int, float] = 60,
        reset: bool = False,
    ):
        """
        :param name: A unique name for this scheduler.
        :param max_calcjobs: Maximum number of calcjobs to run concurrently
        :param max_workflows: Maximum number of top-level workflows to run concurrently
        :param max_processes: Maximum number of processes to run concurrently
        :param poll_interval: Interval in seconds for the fallback polling mechanism
                              (if an AiiDA broadcast is missed).
        """
        self.name = name
        self._max_calcjobs = max_calcjobs
        self._max_workflows = max_workflows
        self._max_processes = max_processes

        self._node: SchedulerNode = None
        self._poll_interval = poll_interval

        # Set up AiiDA profile and event loop policy
        manager = get_manager()
        self._profile = manager.get_profile()
        self._prefix = f"aiida-{self._profile.uuid}"
        set_event_loop_policy()

        # We create our dedicated event loop
        self._loop = asyncio.new_event_loop()

        if reset:
            self.reset()

        # Set up the RabbitMQ communicator and attach the task receiver callback
        self.communicator.add_task_subscriber(self.task_receiver)
        try:
            self.communicator.add_rpc_subscriber(
                self.message_receive, identifier=str(self.node.pk)
            )
            # self.add_cleanup(functools.partial(self.communicator.remove_rpc_subscriber, identifier))
        except kiwipy.TimeoutError:
            raise exceptions.InternalError(
                "Failed to add RPC subscriber to the communicator. "
                "Is the RabbitMQ server running?"
            )

        self.process_controller = get_manager().get_process_controller()

    @property
    def node(self) -> SchedulerNode:
        """Return (and lazily create) the SchedulerNode that stores all data for this scheduler."""
        if self._node is None:
            qb = QueryBuilder()
            qb.append(SchedulerNode, filters={"attributes.name": self.name})
            if qb.count() == 0:
                # Create a new SchedulerNode
                self._node = SchedulerNode(
                    name=self.name,
                    max_calcjobs=self._max_calcjobs,
                    max_workflows=self._max_workflows,
                    max_processes=self._max_processes,
                )
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
                # Check for any stale processes
                self._cleanup_stale_processes()

        return self._node

    def reset(self) -> None:
        """
        Reset the scheduler by clearing the waiting and running process lists.
        """
        self.node.reset()
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
    def queue_name(self):
        return f"{self._prefix}-{self.name}"

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
                task_queue=f"{self._prefix}-{self.name}",
                task_prefetch_count=10000,
                # This is needed because the verdi commands will call this function and when called in unit tests the
                # testing_mode cannot be set.
                testing_mode=self._profile.is_test_profile,
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
            self._add_process(pid)

    def _add_process(self, pk: int | Node) -> None:
        """Add a process to the scheduler."""
        try:
            node = load_node(pk) if not isinstance(pk, Node) else pk
            pk = node.pk
            # Check if the process is already terminated
            if node.is_terminated:
                LOGGER.warning(
                    "Try to add process pk=%d, but it is already terminated.",
                    pk,
                )
                return
            elif node.process_state == ProcessState.CREATED:
                # If the process is in CREATED state, we can continue it
                self.node.append_waiting_process(pk)
            elif node.process_state in [ProcessState.RUNNING, ProcessState.WAITING]:
                # If the process is in RUNNING or WAITING state, we can continue it
                self.node.append_running_process(pk)
                self.call_on_process_finish(pk)
            # The process may already has SCHEDULER_KEY and priority set
            # but we set it again to be sure
            node.base.extras.set(SCHEDULER_KEY, self.name)
            priority = self._compute_priority_for_new_process(pk)
            node.base.extras.set(SCHEDULER_PRIORITY_KEY, priority)
            # Trigger consumption (start more processes if possible)
            self.consume_process_queue()
        except exceptions.NotExistent:
            LOGGER.warning(
                "Received CONTINUE_TASK for pk=%d, but it doesn't exist anymore.",
                pk,
            )

    def message_receive(self, _comm: kiwipy.Communicator, msg: Dict[str, Any]) -> Any:
        """
        Coroutine called when the scheduler receives a message from the communicator

        :param _comm: the communicator that sent the message
        :param msg: the message
        :return: the outcome of processing the message, the return value will be sent back as a response to the sender
        """
        LOGGER.info(
            "Process<%s>: received RPC message: %r",
            self.node.pk,
            msg,
        )

        intent = msg[INTENT_KEY]

        if intent.lower() == Intent.STOP:
            self._loop.call_soon(self.close)

        if intent.lower() == Intent.CONTINUE:
            self._loop.call_soon(self.consume_process_queue)

        if intent.lower() == Intent.PLAY:
            pks = msg[MESSAGE_TEXT_KEY]
            for pk in pks:
                node = self._validate_before_continue(pk)
                if node is not False:
                    self._loop.call_soon(self.continue_process, pk)

        if intent.lower() == Intent.REPAIR:
            pks = msg[MESSAGE_TEXT_KEY]
            self._loop.call_soon(self._repair_processes, pks)

        if intent.lower() == Intent.SET:
            data = msg[MESSAGE_TEXT_KEY]
            identifier = data.get("identifier")
            value = data.get("value")
            if identifier.lower() == "max_calcjobs":
                self.node.max_calcjobs = value
                self.consume_process_queue()
                return self.node.max_calcjobs
            elif identifier.lower() == "max_workflows":
                self.node.max_workflows = value
                self.consume_process_queue()
                return self.node.max_workflows
            elif identifier.lower() == "max_processes":
                self.node.max_processes = value
                self.consume_process_queue()
                return self.node.max_processes
            elif identifier.lower() == "priority":
                processes = data.get("processes")
                for pk in processes:
                    child_node = load_node(pk)
                    child_node.base.extras.set(SCHEDULER_PRIORITY_KEY, value)
            else:
                raise RuntimeError(f"Unknown identifier: {identifier}")

        if intent == Intent.STATUS:
            status_info = {
                "priority": self.node.get_process_priority(),
                "running_process": self.node.running_process,
                "running_workflow": self.node.running_workflow,
                "running_calcjob": self.node.running_calcjob,
            }
            return status_info

        raise RuntimeError("Unknown intent")

    def _compute_priority_for_new_process(self, pid: int | Node) -> int:
        """"""
        from aiida.common.links import LinkType

        # find the parent process
        if isinstance(pid, int):
            node = load_node(pid)
        else:
            node = pid
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

    def _has_capacity(self, node: Node) -> bool:
        """Return True if we have a free slot to run *node* (CalcJob or WorkChain)."""
        if isinstance(node, CalcJobNode):
            return len(self.node.running_calcjob) < self.node.max_calcjobs
        if SchedulerNode.is_top_level_workflow(node):
            return len(self.node.running_workflow) < self.node.max_workflows
        return len(self.node.running_process) < self.node.max_processes

    def _next_waiting_pk(self) -> Optional[int]:
        """
        Return the PK of the waiting process with the **highest** priority
        (largest integer). None if the waiting list is empty.
        """
        priorities = self.node.get_process_priority()
        return max(priorities, key=priorities.get) if priorities else None

    def _log_summary(self):
        LOGGER.info(
            "waiting=%d  process=%d/%d  calcjob=%d/%d, workflow=%d/%d",
            len(self.node.waiting_process),
            len(self.node.running_process),
            self.node.max_processes,
            len(self.node.running_calcjob),
            self.node.max_calcjobs,
            len(self.node.running_workflow),
            self.node.max_workflows,
        )

    def consume_process_queue(self) -> None:
        """
        Try to start as many queued processes as the capacity limits allow.
        Never raises – all errors are logged and the loop continues.
        """

        self._log_summary()

        while True:
            pk = self._next_waiting_pk()
            if pk is None:
                LOGGER.debug("No processes waiting; scheduler idle.")
                break

            node = self._validate_before_continue(pk)
            if node is False:
                # process is already terminated or doesn't exist anymore
                self.node.remove_waiting_process(pk)
                LOGGER.warning(
                    "Try to launch process pk=%d, but it is already terminated or doesn't exist anymore.",
                    pk,
                )
                continue

            if not self._has_capacity(node):
                # queue is full – but there *is* work waiting
                LOGGER.debug("Capacity full; will continue when a slot frees up.")
                break

            try:
                self.continue_process(node)
                self._log_summary()  # live feedback after each launch
            except Exception:  # noqa: BLE001 – keep scheduler alive
                LOGGER.exception("Failed launching pk=%d – skipped.", pk)
                # failed launches are *not* removed from waiting so the user can retry
                break

    def _pop_highest_priority_waiting_process(self) -> Optional[int]:
        """Return the waiting process PK with the largest 'priority' extras, or None if empty."""
        priorites = self.node.get_process_priority()
        # find the process with the highest priority
        if priorites:
            best_pk = max(priorites, key=priorites.get)
            LOGGER.debug(f"Best waiting process pk={best_pk}")
            return best_pk
        else:
            return None

    def _get_priority_for_pid(self, pk: int) -> int:
        """Return the integer priority from the node extras, or a negative number if not set."""
        child_node = load_node(pk)
        return child_node.base.extras.get(SCHEDULER_PRIORITY_KEY, -9999999)

    def _cleanup_stale_processes(self) -> None:
        """
        Clean up processes that are incorrectly listed as running or waiting.

        This function scans the 'running_process' and 'waiting_process' lists for process nodes that are
        no longer active (i.e., terminated or deleted). Terminated processes are removed from the list,
        and callbacks are reattached to still-active running processes to ensure they are monitored properly.

        This is especially useful after scheduler restarts or interruptions where the state might be outdated.
        """

        running_process = self.node.running_process
        if running_process:
            existing_processes = query_existing_processes(running_process)
            non_existing_processes = set(running_process) - set(existing_processes)
            terminated_pks = query_terminated_processes(running_process)
            if terminated_pks:
                LOGGER.info(
                    f"Found {len(terminated_pks)} processes in 'running' "
                    "but it's already terminated. Removing it.",
                )
            if non_existing_processes:
                LOGGER.warning(
                    f"Found {len(non_existing_processes)} processes in 'running' "
                    "but it doesn't exist anymore. Removing it.",
                )
            processes_to_remove = list(non_existing_processes) + list(terminated_pks)
            self.node.remove_running_process(processes_to_remove)
            running_process = set(running_process) - set(processes_to_remove)
            # If we have any processes that are still running, re-attach the callback
            if running_process:
                LOGGER.info(
                    f"Found {len(running_process)} processes in 'running' "
                    "but not terminated. Re-attaching callback.",
                )
                for pk in running_process:
                    self.call_on_process_finish(pk)
        # -----------------------------------------------------------
        waiting_process = self.node.waiting_process
        if waiting_process:
            existing_processes = query_existing_processes(waiting_process)
            non_existing_processes = set(waiting_process) - set(existing_processes)
            terminated_pks = query_terminated_processes(waiting_process)
            if terminated_pks:
                LOGGER.info(
                    f"Found {len(terminated_pks)} processes in 'waiting' "
                    "but it's already terminated. Removing it.",
                )
            if non_existing_processes:
                LOGGER.warning(
                    f"Found {len(non_existing_processes)} processes in 'waiting' "
                    "but it doesn't exist anymore. Removing it.",
                )
            self.node.remove_waiting_process(
                list(non_existing_processes) + list(terminated_pks)
            )

    def _validate_before_continue(self, pk: int) -> bool | Node:
        """
        Validate if the process can be continued.
        This is a placeholder for any validation logic that might be needed.
        """
        # Check if the process is already terminated
        try:
            node = load_node(pk)
            if node.is_terminated:
                LOGGER.warning(
                    "Try to continue process pk=%d, but it is already terminated.",
                    pk,
                )
                return False
            return node
        except exceptions.NotExistent:
            LOGGER.warning(
                "Received CONTINUE_TASK for pk=%d, but it doesn't exist anymore.",
                pk,
            )
            return False

    def continue_process(self, pk: int | Node) -> None:
        """
        Actually continue an AiiDA process in the daemon. Attach a callback
        so that when it finishes, we free up a slot and can launch a new one.
        """
        pk = pk.pk if isinstance(pk, Node) else pk
        try:
            LOGGER.info("Continuing process pk=%d...", pk)
            self.call_on_process_finish(pk)
            # Do not wait for the future's result to avoid blocking
            self.process_controller.continue_process(pk, nowait=False, no_reply=True)
            self.node.append_running_process(pk)
            self.node.remove_waiting_process(pk)
        except Exception:
            LOGGER.exception(
                f"Failed to continue process {pk}. It will be removed from the scheduler.",
            )
            self.node.remove_waiting_process(pk)
            self.node.remove_running_process(pk)

    def call_on_process_finish(self, pk: int) -> None:
        """
        Attach both a broadcast-based subscriber and a fallback polling
        so that when process <pk> terminates, we call `on_process_finished_callback`.
        """
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

    def on_process_finished_callback(self, pk: int) -> None:
        """
        When a process finishes, remove it from the running list, and try
        to consume another from the waiting queue.
        """
        LOGGER.info("[📢] Process pk=%d finished", pk)
        self.node.remove_running_process(pk)
        # Attempt to start more processes if any remain
        self.consume_process_queue()

    def _poll_process(self):
        """Fallback on polling to check if processes are still running.
        This is useful if the RabbitMQ broadcast was missed or not received.
        """
        running_process = self.node.running_process
        if running_process:
            terminated_pks = query_terminated_processes(running_process)
            for pk in terminated_pks:
                LOGGER.info(
                    f"pool_process: {pk} is already terminated (confirmed by polling).",
                )
                self._loop.call_soon(self.on_process_finished_callback, pk)
        self._loop.call_later(self._poll_interval, self._poll_process)

    def start(self) -> None:
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
        # Schedule the first poll
        self._loop.call_later(self._poll_interval, self._poll_process)

        try:
            self.node.is_running = True
            self._loop.run_forever()
        except KeyboardInterrupt:
            LOGGER.info("Scheduler '%s' interrupted by user (Ctrl+C).", self.name)
        finally:
            self._loop.run_until_complete(self._loop.shutdown_asyncgens())
            self._loop.close()
            LOGGER.info("Scheduler '%s' event loop closed.", self.name)

    def close(self) -> None:
        """Close the runner by stopping the loop."""
        self._loop.stop()
        if not self._loop.is_running():
            self._loop.close()
        reset_event_loop_policy()
        self.node.is_running = False
        LOGGER.info("Scheduler '%s' closed.", self.name)

    @classmethod
    def get_scheduler_node(cls, name: str) -> SchedulerNode:
        """
        Return the scheduler node associated with the given name.
        """
        qb = QueryBuilder()
        qb.append(SchedulerNode, filters={"attributes.name": name})
        if qb.count() == 0:
            raise exceptions.NotExistent(f"Scheduler '{name}' does not exist.")
        return qb.first()[0]

    @classmethod
    def get_status(cls, name: str) -> Dict[str, Any]:
        """
        Return the status of the scheduler.
        """
        controller = get_manager().get_process_controller()

        scheduler = cls.get_scheduler_node(name)

        try:
            status = controller._communicator.rpc_send(
                scheduler.pk, {"intent": "status"}
            )
            result = status.result().result()
            return result
        except kiwipy.exceptions.UnroutableError:
            return None

    @classmethod
    def stop_scheduler(cls, name: str) -> None:
        """
        Stop the scheduler with the given name.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        controller._communicator.rpc_send(scheduler.pk, {"intent": "stop"})

    @classmethod
    def set_max_calcjobs(cls, name: str, max_calcjobs: int) -> None:
        """
        Set the maximum number of calcjobs for the scheduler with the given name.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        result = controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "set",
                "message": {"identifier": "max_calcjobs", "value": max_calcjobs},
            },
        )
        return result

    @classmethod
    def set_max_workflows(cls, name: str, max_workflows: int) -> None:
        """
        Set the maximum number of workflows for the scheduler with the given name.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        result = controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "set",
                "message": {"identifier": "max_workflows", "value": max_workflows},
            },
        )
        return result

    @classmethod
    def set_max_processes(cls, name: str, max_processes: int) -> None:
        """
        Set the maximum number of processes for the scheduler with the given name.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        result = controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "set",
                "message": {"identifier": "max_processes", "value": max_processes},
            },
        )
        return result

    @classmethod
    def play_processes(cls, name: str, pks: list, timeout: int = 5) -> None:
        """
        Play processes with the given pks.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "play",
                "message": pks,
            },
        )

    def _repair_processes(self, pks: list) -> None:
        """
        Repair processes with the given pks by continuing them in the scheduler.
        """
        for pk in pks:
            try:
                node = load_node(pk)
                self._add_process(node)
            except exceptions.NotExistent:
                self.node.remove_running_process(pk)
                self.node.remove_waiting_process(pk)
                LOGGER.warning(
                    "Process pk=%d does not exist. Removing it from the scheduler.",
                    pk,
                )

    @classmethod
    def repair_processes(cls, name: str, pks: list, timeout: int = 5) -> None:
        """
        Repair processes with the given pks by continuing them in the scheduler.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "repair",
                "message": pks,
            },
        )

    @classmethod
    def set_process_priority(
        cls, name: str, pks: list, priority: int, timeout: int = 5
    ) -> None:
        """
        Play a process with the given pk.
        """
        scheduler = cls.get_scheduler_node(name)
        controller = get_manager().get_process_controller()
        controller._communicator.rpc_send(
            scheduler.pk,
            {
                "intent": "set",
                "message": {
                    "identifier": "priority",
                    "value": priority,
                    "processes": pks,
                },
            },
        )

    def iterate_tasks(self):
        """Return an iterator over the tasks in the launch queue.
        Note: this will not delete the tasks from the queue, because
        we don't acknowledge them.
        """
        for task in self.communicator._communicator.task_queue(self.queue_name):
            yield task

    def get_process_tasks(self):
        pks = []

        for task in self.iterate_tasks():
            try:
                pks.append(task.body.get("args", {})["pid"])
            except KeyError:
                pass

        return pks
