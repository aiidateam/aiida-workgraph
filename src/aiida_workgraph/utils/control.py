from __future__ import annotations

from aiida.manage import get_manager
from aiida import orm
from aiida.engine.processes import control
from plumpy.process_comms import RemoteProcessThreadController
from typing import Any, Optional, Type


class ControllerWithQueueName(RemoteProcessThreadController):
    def __init__(self, queue_name: str, **kwargs):
        super().__init__(**kwargs)
        self.queue_name = queue_name
        # This is a workaround to patch the communicator to preserve the queue to avoid
        # opening many channels.
        # A better solution is implement it in the manager, so that
        # we can also close the channel when loading another profile.
        if not hasattr(self._communicator, "scheduler_queues"):
            self._communicator.scheduler_queues = {}
        if self.queue_name not in self._communicator.scheduler_queues:
            queue = self._communicator.task_queue(self.queue_name)
            self._communicator.scheduler_queues[self.queue_name] = queue

    @property
    def scheduler_queue(self):
        """Return the scheduler queue."""
        return self._communicator.scheduler_queues[self.queue_name]

    def task_send(self, message: Any, no_reply: bool = False) -> Optional[Any]:
        """
        Send a task to be performed using the communicator

        :param message: the task message
        :param no_reply: if True, this call will be fire-and-forget, i.e. no return value
        :return: the response from the remote side (if no_reply=False)
        """
        result = self.scheduler_queue.task_send(message, no_reply=no_reply)
        return result


def create_task_action(
    pk: int,
    tasks: list,
    action: str = "pause",
    timeout: int = 5,
    wait: bool = False,
):
    """Send task action to Process."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "catalog": "task", "action": action, "tasks": tasks}
    controller._communicator.rpc_send(pk, message)


def continue_process_in_scheduler(
    pk: int | orm.Node,
    scheduler_name: str = "test_scheduler",
):
    """Send workgraph task to scheduler."""

    manager = get_manager()
    profile = manager.get_profile()
    queue_name = f"aiida-{profile.uuid}-{scheduler_name}"
    controller = ControllerWithQueueName(
        queue_name=queue_name, communicator=manager.get_communicator()
    )
    if isinstance(pk, orm.Node):
        node = pk
        pk = node.pk
    else:
        node = orm.load_node(pk)
    controller.continue_process(pk, nowait=False)
    node.base.extras.set("_scheduler", scheduler_name)


def play_process_in_scheduler(scheduler: int, pk: int | orm.Node):
    controller = get_manager().get_process_controller()
    controller._communicator.rpc_send(scheduler, {"intent": "play", "message": pk})


def stop_scheduler(scheduler: int):
    """Stop scheduler."""
    controller = get_manager().get_process_controller()
    controller._communicator.rpc_send(scheduler, {"intent": "stop"})


def get_scheduler_status(scheduler: int):
    """Get scheduler status."""
    controller = get_manager().get_process_controller()
    status = controller._communicator.rpc_send(scheduler, {"intent": "status"})
    result = status.result().result()
    return result


def get_task_runtime_info(node, name: str, key: str) -> str:
    """Get task state info from attributes."""
    from aiida_workgraph.orm.utils import deserialize_safe

    if key == "process":
        value = deserialize_safe(node.task_processes.get(name, ""))
    elif key == "state":
        value = node.task_states.get(name, "")
    elif key == "action":
        value = node.task_actions.get(name, "")
    return value


def pause_tasks(pk: int, tasks: list[str], timeout: int = 5, wait: bool = False):
    """Pause task."""
    node = orm.load_node(pk)
    if node.is_finished:
        message = "WorkGraph is finished. Cannot pause tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            if get_task_runtime_info(node, name, "state") == "PLANNED":
                create_task_action(pk, tasks, action="pause")
            elif get_task_runtime_info(node, name, "state") == "RUNNING":
                try:
                    control.pause_processes(
                        [get_task_runtime_info(node, name, "process")],
                        all_entries=None,
                        timeout=5,
                        wait=False,
                    )
                except Exception as e:
                    print(f"Pause task {name} failed: {e}")
    return True, ""


def play_tasks(pk: int, tasks: list, timeout: int = 5, wait: bool = False):
    node = orm.load_node(pk)
    if node.is_finished:
        message = "WorkGraph is finished. Cannot kill tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            if get_task_runtime_info(node, name, "state") == "PLANNED":
                create_task_action(pk, tasks, action="play")
            elif get_task_runtime_info(node, name, "state") in ["CREATED", "PAUSED"]:
                try:
                    control.play_processes(
                        [get_task_runtime_info(node, name, "process")],
                        all_entries=None,
                        timeout=5,
                        wait=False,
                    )
                except Exception as e:
                    print(f"Play task {name} failed: {e}")
            elif get_task_runtime_info(node, name, "process").is_finished:
                raise ValueError(f"Task {name} is already finished.")
    return True, ""


def kill_tasks(pk: int, tasks: list, timeout: int = 5, wait: bool = False):
    node = orm.load_node(pk)
    if node.is_finished:
        message = "WorkGraph is finished. Cannot kill tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        print("tasks", tasks)
        for name in tasks:
            state = get_task_runtime_info(node, name, "state")
            process = get_task_runtime_info(node, name, "process")
            print("state", state)
            print("process", process)
            if state == "PLANNED":
                create_task_action(pk, tasks, action="skip")
            elif state in [
                "CREATED",
                "RUNNING",
                "WAITING",
                "PAUSED",
            ]:
                if get_task_runtime_info(node, name, "process") is None:
                    print(f"Task {name} is not a AiiDA process.")
                    create_task_action(pk, tasks, action="kill")
                else:
                    try:
                        control.kill_processes(
                            [get_task_runtime_info(node, name, "process")],
                            all_entries=None,
                            timeout=5,
                            wait=False,
                        )
                    except Exception as e:
                        print(f"Kill task {name} failed: {e}")
    return True, ""


def reset_tasks(pk: int, tasks: list) -> None:
    """Reset tasks
    Args:
        tasks (list): a list of task names.
    """
    node = orm.load_node(pk)
    if node.is_finished:
        message = "WorkGraph is finished. Cannot kill tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            create_task_action(pk, tasks, action="reset")

    return True, ""


def submit_to_scheduler(
    process_class, inputs: dict, scheduler: str | None = None
) -> orm.Node:
    """Submit a process to the scheduler."""
    from aiida.manage import manager
    from aiida.engine.utils import instantiate_process
    from aiida.engine import submit
    from aiida_workgraph.utils.control import continue_process_in_scheduler

    if scheduler is not None:
        runner = manager.get_manager().get_runner()
        process_inited = instantiate_process(runner, process_class, **inputs)
        process_inited.close()
        process_inited.runner.persister.save_checkpoint(process_inited)
        node = process_inited.node
        continue_process_in_scheduler(node.pk, scheduler_name=scheduler)
    else:
        node = submit(process_class, **inputs)

    return node


def submit_to_scheduler_inside_workchain(
    self,
    process: Type["Process"],
    inputs: dict[str, Any] | None = None,
    **kwargs,
) -> orm.ProcessNode:
    """Submit a process inside the workchain to the scheduler.
    :param process: The process class.
    :param inputs: The dictionary of process inputs.
    :return: The process node.
    """
    from aiida.engine import utils
    from aiida_workgraph.utils.control import continue_process_in_scheduler

    scheduler = self.node.base.extras.get("_scheduler", None)
    if scheduler:
        inputs = utils.prepare_inputs(inputs, **kwargs)
        process_inited = self.runner.instantiate_process(process, **inputs)
        self.runner.persister.save_checkpoint(process_inited)
        process_inited.close()
        continue_process_in_scheduler(process_inited.node, scheduler)
        return process_inited.node
    else:
        return self.runner.submit(process, inputs, **kwargs)
