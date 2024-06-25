from aiida.manage import get_manager
from aiida.orm import ProcessNode


def create_task_action(
    process: ProcessNode,
    tasks: list,
    action: str = "pause",
    timeout: int = 5,
    wait: bool = False,
):
    """Send task action to Process."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "catalog": "task", "action": action, "tasks": tasks}
    controller._communicator.rpc_send(process.pk, message)
