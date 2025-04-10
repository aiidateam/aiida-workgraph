from aiida.orm import Data
from aiida.common.lang import classproperty
from aiida.orm.fields import add_field
from aiida.orm.utils.mixins import Sealable
from typing import Optional, Tuple


class SchedulerNode(Sealable, Data):
    """"""

    WAITING_PROCESS = "waiting_process"
    RUNNING_PROCESS = "running_process"
    FINISHED_PROCESS = "finished_process"
    MAXIMUM_RUNNING_PROCESSES = "maximum_running_processes"

    __qb_fields__ = [
        add_field(
            "name",
            dtype=Optional[str],
            doc="Name of the scheduler node",
        ),
        add_field(
            WAITING_PROCESS,
            dtype=Optional[list],
            doc="List of waiting processes",
        ),
        add_field(
            RUNNING_PROCESS,
            dtype=Optional[list],
            doc="List of running processes",
        ),
        add_field(
            FINISHED_PROCESS,
            dtype=Optional[list],
            doc="List of finished processes",
        ),
        add_field(
            MAXIMUM_RUNNING_PROCESSES,
            dtype=Optional[int],
            doc="Maximum number of running processes",
        ),
    ]

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # noqa: N805
        return super()._updatable_attributes + (
            cls.WAITING_PROCESS,
            cls.RUNNING_PROCESS,
            cls.FINISHED_PROCESS,
            cls.MAXIMUM_RUNNING_PROCESSES,
        )

    @property
    def name(self) -> str:
        """Return the name of the scheduler node."""
        return self.base.attributes.get("name", "scheduler")

    @name.setter
    def name(self, value: str) -> None:
        """Set the name of the scheduler node."""
        self.base.attributes.set("name", value)

    @property
    def waiting_process(self) -> list:
        """Return the list of waiting processes."""
        return self.base.attributes.get(self.WAITING_PROCESS, [])

    @waiting_process.setter
    def waiting_process(self, value: list) -> None:
        """Set the list of waiting processes."""
        self.base.attributes.set(self.WAITING_PROCESS, value)

    def append_waiting_process(self, pk: int) -> None:
        waiting_process = self.base.attributes.get(self.WAITING_PROCESS, [])
        if pk not in waiting_process:
            waiting_process.append(pk)
            self.waiting_process = waiting_process

    def pop_waiting_process(self) -> None:
        waiting_process = self.base.attributes.get(self.WAITING_PROCESS, [])
        if waiting_process:
            pk = waiting_process.pop(0)
            self.waiting_process = waiting_process
            return pk
        return None

    @property
    def running_processes(self) -> list:
        """Return the list of running processes."""
        return self.base.attributes.get(self.RUNNING_PROCESS, [])

    @running_processes.setter
    def running_processes(self, value: list) -> None:
        """Set the list of running processes."""
        self.base.attributes.set(self.RUNNING_PROCESS, value)

    def append_running_process(self, pk: int) -> None:
        running_process = self.base.attributes.get(self.RUNNING_PROCESS, [])
        if pk not in running_process:
            running_process.append(pk)
            self.running_processes = running_process

    def remove_running_process(self, pk: int) -> None:
        running_process = self.base.attributes.get(self.RUNNING_PROCESS, [])
        if pk in running_process:
            running_process.remove(pk)
            self.running_processes = running_process

    @property
    def maxium_running_processes(self) -> int:
        """Return the maximum number of running processes."""
        return self.base.attributes.get(self.MAXIMUM_RUNNING_PROCESSES, 100000)

    @maxium_running_processes.setter
    def maxium_running_processes(self, value: int) -> None:
        """Set the maximum number of running processes."""
        self.base.attributes.set(self.MAXIMUM_RUNNING_PROCESSES, value)
