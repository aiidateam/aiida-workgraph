"""Module with `Node` sub class for work processes."""
from typing import Optional, Tuple

from aiida.common.lang import classproperty

from aiida.orm.nodes.process.workflow.workchain import WorkChainNode

__all__ = ("WorkGraphNode",)


class WorkGraphNode(WorkChainNode):
    """ORM class for all nodes representing the execution of a WorkGraph."""

    TASK_STATES_KEY = "task_states"
    TASK_PROCESSES_KEY = "task_processes"
    TASK_ACTIONS_KEY = "task_actions"
    WORKGRAPH_DATA_KEY = "workgraph_data"

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # type: ignore
        # pylint: disable=no-self-argument
        return super()._updatable_attributes + (
            cls.TASK_STATES_KEY,
            cls.TASK_PROCESSES_KEY,
            cls.TASK_ACTIONS_KEY,
        )

    @property
    def task_states(self) -> Optional[str]:
        """
        Return the task state info

        :returns: string representation of the task state info
        """
        return self.base.attributes.get(self.TASK_STATES_KEY, {})

    def set_task_states(self, task_states: str) -> None:
        """
        Set the task state info

        :param state: string representation of the task state info
        """
        return self.base.attributes.set(self.TASK_STATES_KEY, task_states)

    def get_task_state(self, task_name: str) -> Optional[str]:
        """
        Return the task state info

        :returns: string representation of the task state info
        """
        return self.task_states.get(task_name, "")

    def set_task_state(self, task_name: str, task_state: str) -> Optional[str]:
        """
        Set the task state info

        :returns: string representation of the task state info
        """
        task_states = self.task_states
        task_states[task_name] = task_state
        return self.set_task_states(task_states)

    @property
    def task_processes(self) -> Optional[str]:
        """
        Return the task process info

        :returns: string representation of the task process info
        """
        return self.base.attributes.get(self.TASK_PROCESSES_KEY, {})

    def set_task_processes(self, task_processes: str) -> None:
        """
        Set the task process info

        :param process: string representation of the task process info
        """
        return self.base.attributes.set(self.TASK_PROCESSES_KEY, task_processes)

    def get_task_process(self, task_name: str) -> Optional[str]:
        """
        Return the task state info

        :returns: string representation of the task state info
        """
        return self.task_processes.get(task_name, None)

    def set_task_process(self, task_name: str, task_process: str) -> Optional[str]:
        """
        Set the task state info

        :returns: string representation of the task state info
        """
        task_processes = self.task_processes
        task_processes[task_name] = task_process
        return self.set_task_processes(task_processes)

    @property
    def task_actions(self) -> Optional[str]:
        """
        Return the task action info

        :returns: string representation of the task action info
        """
        return self.base.attributes.get(self.TASK_ACTIONS_KEY, {})

    def set_task_actions(self, task_actions: str) -> None:
        """
        Set the task action info

        :param action: string representation of the task action info
        """
        return self.base.attributes.set(self.TASK_ACTIONS_KEY, task_actions)

    def get_task_action(self, task_name: str) -> Optional[str]:
        """
        Return the task action info

        :returns: string representation of the task action info
        """
        return self.task_actions.get(task_name, "")

    def set_task_action(self, task_name: str, task_action: str) -> Optional[str]:
        """
        Set the task action info

        :returns: string representation of the task action info
        """
        task_actions = self.task_actions
        task_actions[task_name] = task_action
        return self.set_task_actions(task_actions)

    @property
    def workgraph_data(self) -> Optional[str]:
        """
        Return the workgraph data

        :returns: string representation of the workgraph data
        """
        return self.base.attributes.get(self.WORKGRAPH_DATA_KEY, None)

    def set_workgraph_data(self, workgraph_data: dict) -> None:
        """
        Set the workgraph data

        :param data: string representation of the workgraph data
        """
        return self.base.attributes.set(self.WORKGRAPH_DATA_KEY, workgraph_data)
