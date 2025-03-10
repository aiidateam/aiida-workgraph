from __future__ import annotations
from typing import Dict


class TaskActionManager:
    """
    Handles externally-triggered task actions (RESET, PAUSE, PLAY, SKIP, KILL)
    to change task states or runtime behavior.
    """

    def __init__(self, state_manager, logger, process):
        """
        :param state_manager: A reference to TaskStateManager for updating states.
        :param logger: A logger instance.
        """
        self.state_manager = state_manager
        self.logger = logger
        self.process = process

    def apply_task_actions(self, msg: Dict) -> None:
        """
        Apply task actions to the workgraph based on user or external messages.
        :param msg: { "action": <str>, "tasks": <List[str]> }
        """
        action = msg["action"].upper()
        tasks = msg["tasks"]
        self.process.report(f"Action: {action}. Tasks: {tasks}")

        if action == "RESET":
            for name in tasks:
                self.state_manager.reset_task(name)
        elif action == "PAUSE":
            for name in tasks:
                self.pause_task(name)
        elif action == "PLAY":
            for name in tasks:
                self.play_task(name)
        elif action == "SKIP":
            for name in tasks:
                self.skip_task(name)
        elif action == "KILL":
            for name in tasks:
                self.kill_task(name)

    def pause_task(self, name: str) -> None:
        """Mark the task to be paused."""
        self.state_manager.set_task_runtime_info(name, "action", "PAUSE")
        self.process.report(f"Task {name} action: PAUSE.")

    def play_task(self, name: str) -> None:
        """Remove the PAUSE action on a task."""
        self.state_manager.set_task_runtime_info(name, "action", "")
        self.process.report(f"Task {name} action: PLAY.")

    def skip_task(self, name: str) -> None:
        """Force a task to be SKIPPED."""
        self.state_manager.set_task_runtime_info(name, "state", "SKIPPED")
        self.process.report(f"Task {name} action: SKIP.")

    def kill_task(self, name: str) -> None:
        """
        KILL a running task. Typically used for AWAITABLE or MONITOR tasks
        to cancel the underlying async future.
        """
        state = self.state_manager.get_task_runtime_info(name, "state")
        if state == "RUNNING":
            task = self.process.wg.tasks[name]
            node_type = task.node_type.upper()
            if node_type in ["AWAITABLE", "MONITOR"]:
                awaitable_manager = self.state_manager.awaitable_manager
                awaitable_target = awaitable_manager.not_persisted_awaitables.get(name)
                if awaitable_target:
                    try:
                        awaitable_target.cancel()
                        self.state_manager.set_task_runtime_info(
                            name, "state", "KILLED"
                        )
                        self.process.report(f"Task {name} was KILLED.")
                    except Exception as e:
                        self.logger.error(f"Error in killing task {name}: {e}")
                else:
                    self.logger.warning(f"No active awaitable found for task {name}.")
