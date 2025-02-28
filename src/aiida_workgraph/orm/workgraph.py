"""Module with `Node` sub class for work processes."""
from typing import Optional, Tuple
import logging
from aiida.common.lang import classproperty

from aiida.orm.nodes.process.workflow.workchain import WorkChainNode

__all__ = ("WorkGraphNode",)


class WorkGraphNode(WorkChainNode):
    """ORM class for all nodes representing the execution of a WorkGraph."""

    TASK_STATES_KEY = "task_states"
    TASK_PROCESSES_KEY = "task_processes"
    TASK_ACTIONS_KEY = "task_actions"
    TASK_EXECUTORS_KEY = "task_executors"
    TASK_ERROR_HANDLERS_KEY = "task_error_handlers"
    WORKGRAPH_DATA_KEY = "workgraph_data"
    WORKGRAPH_DATA_SHORT_KEY = "workgraph_data_short"
    WORKGRAPH_ERROR_HANDLERS_KEY = "workgraph_error_handlers"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Use the same logger as WorkChainNode, this ensures the log level is set correctly
        self._logger = logging.getLogger(
            "aiida.orm.nodes.process.workflow.workchain.WorkChainNode"
        )

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # type: ignore
        # pylint: disable=no-self-argument
        return super()._updatable_attributes + (
            cls.WORKGRAPH_DATA_KEY,
            cls.WORKGRAPH_DATA_SHORT_KEY,
            cls.WORKGRAPH_ERROR_HANDLERS_KEY,
            cls.TASK_STATES_KEY,
            cls.TASK_PROCESSES_KEY,
            cls.TASK_ACTIONS_KEY,
            cls.TASK_EXECUTORS_KEY,
            cls.TASK_ERROR_HANDLERS_KEY,
        )

    @property
    def task_states(self) -> Optional[dict]:
        """
        Return the task state info

        :returns: dict of all task states info
        """
        return self.base.attributes.get(self.TASK_STATES_KEY, {})

    @task_states.setter
    def task_states(self, task_states: dict) -> None:
        """
        Set the task state info

        :param state: dict of all task states
        """
        self.base.attributes.set(self.TASK_STATES_KEY, task_states)

    def get_task_state(self, task_name: str) -> Optional[str]:
        """
        Return the task state info

        :returns: the task state info
        """
        return self.task_states.get(task_name, "")

    def set_task_state(self, task_name: str, task_state: str) -> None:
        """Set the task state info"""

        task_states = self.task_states
        task_states[task_name] = task_state
        self.task_states = task_states

    @property
    def task_processes(self) -> Optional[dict]:
        """
        Return the task process info

        :returns: dict of all task process info
        """
        return self.base.attributes.get(self.TASK_PROCESSES_KEY, {})

    @task_processes.setter
    def task_processes(self, task_processes: dict) -> None:
        """
        Set the task process info

        :param task_processes: dict representation of the task process info
        """
        return self.base.attributes.set(self.TASK_PROCESSES_KEY, task_processes)

    def get_task_process(self, task_name: str) -> Optional[str]:
        """
        Return the task process info

        :returns: the task process info
        """
        return self.task_processes.get(task_name, None)

    def set_task_process(self, task_name: str, task_process: str) -> None:
        """Set the task state info"""

        task_processes = self.task_processes
        task_processes[task_name] = task_process
        self.task_processes = task_processes

    @property
    def task_actions(self) -> Optional[dict]:
        """
        Return the task action info

        :returns: string representation of the task action info
        """
        return self.base.attributes.get(self.TASK_ACTIONS_KEY, {})

    @task_actions.setter
    def task_actions(self, task_actions: dict) -> None:
        """
        Set the task action info

        :param task_actions: dict representation of the task action info
        """
        self.base.attributes.set(self.TASK_ACTIONS_KEY, task_actions)

    def get_task_action(self, task_name: str) -> Optional[str]:
        """
        Return the task action info

        :returns: the task action info
        """
        return self.task_actions.get(task_name, "")

    def set_task_action(self, task_name: str, task_action: str) -> None:
        """Set the task action info"""
        task_actions = self.task_actions
        task_actions[task_name] = task_action
        self.task_actions = task_actions

    @property
    def task_executors(self) -> Optional[dict]:
        """
        Return the task executors

        :returns: dict representation of the task executors
        """
        return self.base.attributes.get(self.TASK_EXECUTORS_KEY, {})

    @task_executors.setter
    def task_executors(self, task_executors: dict) -> None:
        """
        Set the task executors

        :param task_executors: dict representation of the task executors
        """
        self.base.attributes.set(self.TASK_EXECUTORS_KEY, task_executors)

    @property
    def task_error_handlers(self) -> Optional[dict]:
        """
        Return the task error handlers

        :returns: dict representation of the task error handlers
        """
        return self.base.attributes.get(self.TASK_ERROR_HANDLERS_KEY, {})

    @task_error_handlers.setter
    def task_error_handlers(self, task_error_handlers: dict) -> None:
        """
        Set the task error handlers

        :param task_error_handlers: dict representation of the task error handlers
        """
        self.base.attributes.set(self.TASK_ERROR_HANDLERS_KEY, task_error_handlers)

    @property
    def workgraph_data(self) -> Optional[dict]:
        """
        Return the workgraph data

        :returns: dict representation of the workgraph data
        """
        return self.base.attributes.get(self.WORKGRAPH_DATA_KEY, None)

    @workgraph_data.setter
    def workgraph_data(self, workgraph_data: dict) -> None:
        """
        Set the workgraph data

        :param workgraph_data: dict representation of the workgraph data
        """
        return self.base.attributes.set(self.WORKGRAPH_DATA_KEY, workgraph_data)

    @property
    def workgraph_data_short(self) -> Optional[dict]:
        """
        Return the short workgraph data, which only has minimal information to be displayed in the GUI.

        :returns: dict representation of the short workgraph data
        """
        return self.base.attributes.get(self.WORKGRAPH_DATA_SHORT_KEY, None)

    @workgraph_data_short.setter
    def workgraph_data_short(self, workgraph_data_short: dict) -> None:
        """
        Set the short workgraph data.

        :param workgraph_data_short: dict representation of the short workgraph data
        """
        return self.base.attributes.set(
            self.WORKGRAPH_DATA_SHORT_KEY, workgraph_data_short
        )

    @property
    def workgraph_error_handlers(self) -> Optional[dict]:
        """
        Return the workgraph error handlers

        :returns: dict representation of the error handlers
        """
        return self.base.attributes.get(self.WORKGRAPH_ERROR_HANDLERS_KEY, None)

    @workgraph_error_handlers.setter
    def workgraph_error_handlers(self, workgraph_error_handlers: dict) -> None:
        """
        Set the workgraph error handlers

        :param workgraph_error_handlers: dict representation of the error handlers
        """
        self.base.attributes.set(
            self.WORKGRAPH_ERROR_HANDLERS_KEY, workgraph_error_handlers
        )
