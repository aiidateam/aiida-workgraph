from aiida.manage import get_manager


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
