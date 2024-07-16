from aiida.manage import get_manager
from aiida import orm
from aiida.engine.processes import control


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


def get_task_state_info(node, name: str, key: str) -> str:
    """Get task state info from base.extras."""
    from aiida.orm.utils.serialize import deserialize_unsafe

    if key == "process":
        value = deserialize_unsafe(node.base.extras.get(f"_task_{key}_{name}", ""))
    else:
        value = node.base.extras.get(f"_task_{key}_{name}", "")
    return value


def set_task_state_info(node, name: str, key: str, value: any) -> None:
    """Set task state info to base.extras."""
    from aiida.orm.utils.serialize import serialize

    if key == "process":
        node.base.extras.set(f"_task_{key}_{name}", serialize(value))
    else:
        node.base.extras.set(f"_task_{key}_{name}", value)


def pause_tasks(pk: int, tasks: list, timeout: int = 5, wait: bool = False):
    """Pause task."""
    node = orm.load_node(pk)
    if node.is_finished:
        message = "Process is finished. Cannot pause tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            if get_task_state_info(node, name, "state") == "PLANNED":
                set_task_state_info(node, name, "action", "pause")
            elif get_task_state_info(node, name, "state") == "RUNNING":
                try:
                    control.pause_processes(
                        [get_task_state_info(node, name, "process")],
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
        message = "Process is finished. Cannot pause tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            if get_task_state_info(node, name, "state") == "PLANNED":
                set_task_state_info(node, name, "action", None)
            elif get_task_state_info(node, name, "state") in ["CREATED", "PAUSED"]:
                try:
                    control.play_processes(
                        [get_task_state_info(node, name, "process")],
                        all_entries=None,
                        timeout=5,
                        wait=False,
                    )
                except Exception as e:
                    print(f"Play task {name} failed: {e}")
            elif get_task_state_info(node, name, "process").is_finished:
                raise ValueError(f"Task {name} is already finished.")
    return True, ""


def kill_tasks(pk: int, tasks: list, timeout: int = 5, wait: bool = False):
    node = orm.load_node(pk)
    if node.is_finished:
        message = "Process is finished. Cannot pause tasks."
        print(message)
        return False, message
    elif node.process_state.value.upper() in [
        "CREATED",
        "RUNNING",
        "WAITING",
        "PAUSED",
    ]:
        for name in tasks:
            if get_task_state_info(node, name, "state") == "PLANNED":
                set_task_state_info(node, name, "action", "skip")
            elif get_task_state_info(node, name, "state") in [
                "CREATED",
                "RUNNING",
                "WAITING",
                "PAUSED",
            ]:
                try:
                    control.kill_processes(
                        [get_task_state_info(node, name, "process")],
                        all_entries=None,
                        timeout=5,
                        wait=False,
                    )
                except Exception as e:
                    print(f"Kill task {name} failed: {e}")
            elif get_task_state_info(node, name, "process").is_finished:
                raise ValueError(f"Task {name} is already finished.")
    return True, ""
