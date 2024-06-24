from aiida.manage import get_manager
from aiida.orm import ProcessNode


def pause_task(process: ProcessNode, task: str, timeout: int = 5, wait: bool = False):
    """Pause task."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "message": f"task,{task}:pause"}
    controller._communicator.rpc_send(process.pk, message)


def play_task(process: ProcessNode, task: str, timeout: int = 5, wait: bool = False):
    """Play task."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "message": f"task,{task}:play"}
    controller._communicator.rpc_send(process.pk, message)


def skip_task(process: ProcessNode, task: str, timeout: int = 5, wait: bool = False):
    """Skip task."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "message": f"task,{task}:skip"}
    controller._communicator.rpc_send(process.pk, message)


def reset_task(process: ProcessNode, task: str, timeout: int = 5, wait: bool = False):
    """Reset task."""

    controller = get_manager().get_process_controller()
    message = {"intent": "custom", "message": f"task,{task}:reset"}
    controller._communicator.rpc_send(process.pk, message)
