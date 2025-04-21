from aiida.orm import Data, load_node, Node, CalcJobNode, QueryBuilder, ProcessNode
from aiida.common.lang import classproperty
from aiida.orm.fields import add_field
from aiida.orm.utils.mixins import Sealable
from typing import Optional, Tuple, Union


class SchedulerNode(Sealable, Data):
    """"""

    IS_RUNNING = "is_running"
    WAITING_PROCESS = "waiting_process"
    RUNNING_PROCESS = "running_process"
    RUNNING_CALCJOB = "running_calcjob"
    RUNNING_WORKFLOW = "running_workflow"
    FINISHED_PROCESS = "finished_process"
    MAX_PROCESSES = "max_processes"
    MAX_CALCJOBS = "max_calcjobs"
    MAX_WORKFLOWS = "max_workflows"
    NEXT_PRIORITY = "next_priority"

    __qb_fields__ = [
        add_field(
            "name",
            dtype=Optional[str],
            doc="Name of the scheduler node",
        ),
        add_field(
            "is_running",
            dtype=Optional[bool],
            doc="Is the scheduler running",
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
            RUNNING_WORKFLOW,
            dtype=Optional[list],
            doc="List of running workflows",
        ),
        add_field(
            RUNNING_CALCJOB,
            dtype=Optional[list],
            doc="List of running calcjobs",
        ),
        add_field(
            FINISHED_PROCESS,
            dtype=Optional[list],
            doc="List of finished processes",
        ),
        add_field(
            MAX_PROCESSES,
            dtype=Optional[int],
            doc="Maximum number of processes",
        ),
        add_field(
            MAX_CALCJOBS,
            dtype=Optional[int],
            doc="Maximum number of running calcjobs",
        ),
        add_field(
            MAX_WORKFLOWS,
            dtype=Optional[int],
            doc="Maximum number of workflows",
        ),
        add_field(
            NEXT_PRIORITY,
            dtype=Optional[int],
            doc="Next priority for the process",
        ),
    ]

    def __init__(
        self,
        name: str = "scheduler",
        max_calcjobs: Optional[int] = None,
        max_workflows: Optional[int] = None,
        max_processes: Optional[int] = None,
        **kwargs,
    ):
        """Initialise a ``SchedulerNode`` instance."""
        super().__init__(**kwargs)
        self.name = name
        self.max_calcjobs = 100 if max_calcjobs is None else max_calcjobs
        self.max_workflows = 100 if max_workflows is None else max_workflows
        self.max_processes = 10000 if max_processes is None else max_processes

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # noqa: N805
        return super()._updatable_attributes + (
            cls.IS_RUNNING,
            cls.WAITING_PROCESS,
            cls.RUNNING_PROCESS,
            cls.RUNNING_WORKFLOW,
            cls.RUNNING_CALCJOB,
            cls.FINISHED_PROCESS,
            cls.MAX_PROCESSES,
            cls.MAX_CALCJOBS,
            cls.MAX_WORKFLOWS,
            cls.NEXT_PRIORITY,
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
    def is_running(self) -> bool:
        """Return if the scheduler is running."""
        return self.base.attributes.get("is_running", False)

    @is_running.setter
    def is_running(self, value: bool) -> None:
        """Set if the scheduler is running."""
        self.base.attributes.set("is_running", value)

    @property
    def waiting_process(self) -> list:
        """Return the list of waiting processes."""
        return self.base.attributes.get(self.WAITING_PROCESS, [])

    @waiting_process.setter
    def waiting_process(self, value: list) -> None:
        """Set the list of waiting processes."""
        self.base.attributes.set(self.WAITING_PROCESS, value)

    def append_waiting_process(self, pk: int) -> None:
        waiting_process = self.waiting_process
        if pk not in waiting_process:
            waiting_process.append(pk)
            self.waiting_process = waiting_process

    def remove_waiting_process(self, pks: Union[int, list]) -> None:
        if isinstance(pks, int):
            pks = [pks]

        waiting_process = set(self.waiting_process)
        waiting_process.difference_update(pks)
        self.waiting_process = list(waiting_process)

    @property
    def running_process(self) -> list:
        """Return the list of running processes."""
        return self.base.attributes.get(self.RUNNING_PROCESS, [])

    @running_process.setter
    def running_process(self, value: list) -> None:
        """Set the list of running processes."""
        self.base.attributes.set(self.RUNNING_PROCESS, value)

    @property
    def running_workflow(self) -> list:
        """Return the list of running workflows."""
        return self.base.attributes.get(self.RUNNING_WORKFLOW, [])

    @running_workflow.setter
    def running_workflow(self, value: list) -> None:
        """Set the list of running workflows."""
        self.base.attributes.set(self.RUNNING_WORKFLOW, value)

    def append_running_workflow(self, pk: int) -> None:
        running_workflow = self.running_workflow
        if pk not in running_workflow:
            running_workflow.append(pk)
            self.running_workflow = running_workflow

    def remove_running_workflow(self, pks: Union[int, list]) -> None:
        if isinstance(pks, int):
            pks = [pks]

        running_workflow = set(self.running_workflow)
        running_workflow.difference_update(pks)
        self.running_workflow = list(running_workflow)

    @property
    def running_calcjob(self) -> list:
        """Return the list of running calcjobs."""
        return self.base.attributes.get(self.RUNNING_CALCJOB, [])

    @running_calcjob.setter
    def running_calcjob(self, value: list) -> None:
        """Set the list of running calcjobs."""
        self.base.attributes.set(self.RUNNING_CALCJOB, value)

    def append_running_calcjob(self, pk: int) -> None:
        running_calcjob = self.running_calcjob
        if pk not in running_calcjob:
            running_calcjob.append(pk)
            self.running_calcjob = running_calcjob

    def remove_running_calcjob(self, pks: Union[int, list]) -> None:
        if isinstance(pks, int):
            pks = [pks]

        running_calcjob = set(self.running_calcjob)
        running_calcjob.difference_update(pks)
        self.running_calcjob = list(running_calcjob)

    @classmethod
    def is_top_level_workflow(cls, node: Node) -> bool:
        """Check if the process is a top level workflow."""
        from aiida.common.links import LinkType

        links = node.base.links.get_incoming()
        # find the parent node
        parent_nodes = [
            link.node
            for link in links
            if link.link_type in [LinkType.CALL_CALC, LinkType.CALL_WORK]
        ]
        return len(parent_nodes) == 0

    def append_running_process(self, pk: int) -> None:

        running_process = self.running_process
        if pk not in running_process:
            running_process.append(pk)
            self.running_process = running_process
            # check if the process is a calcjob
            node = load_node(pk)
            if isinstance(node, CalcJobNode):
                self.append_running_calcjob(pk)
            else:
                if self.is_top_level_workflow(node):
                    self.append_running_workflow(pk)

    def remove_running_process(self, pks: Union[int, list]) -> None:
        if isinstance(pks, int):
            pks = [pks]

        running_process = set(self.running_process)
        running_process.difference_update(pks)
        self.running_process = list(running_process)
        # Also remove from running_calcjob
        self.remove_running_calcjob(pks)
        # Also remove from running_workflow
        self.remove_running_workflow(pks)

    @property
    def max_workflows(self) -> int:
        """Return the maximum number of workflows."""
        return self.base.attributes.get(self.MAX_WORKFLOWS, 100)

    @max_workflows.setter
    def max_workflows(self, value: int) -> None:
        """Set the maximum number of workflows."""
        self.base.attributes.set(self.MAX_WORKFLOWS, value)

    @property
    def max_calcjobs(self) -> int:
        """Return the maximum number of running processes."""
        return self.base.attributes.get(self.MAX_CALCJOBS, 100)

    @max_calcjobs.setter
    def max_calcjobs(self, value: int) -> None:
        """Set the maximum number of running processes."""
        self.base.attributes.set(self.MAX_CALCJOBS, value)

    @property
    def max_processes(self) -> int:
        """Return the maximum number of processes."""
        return self.base.attributes.get(self.MAX_PROCESSES, 10000)

    @max_processes.setter
    def max_processes(self, value: int) -> None:
        """Set the maximum number of processes."""
        self.base.attributes.set(self.MAX_PROCESSES, value)

    @property
    def next_priority(self) -> int:
        """Return the next priority for the process."""
        return self.base.attributes.get(self.NEXT_PRIORITY, 0)

    @next_priority.setter
    def next_priority(self, value: int) -> None:
        """Set the next priority for the process."""
        self.base.attributes.set(self.NEXT_PRIORITY, value)

    def get_process_priority(self) -> dict:
        waiting_list = self.waiting_process
        if not waiting_list:
            return {}

        qb = QueryBuilder()
        qb.append(
            ProcessNode,
            filters={"id": {"in": waiting_list}},
            project=["id", "extras._scheduler_priority"],
        )
        results = qb.all()
        priorites = {result[0]: result[1] for result in results}
        # For a unknown reason, the query returns None for some of the nodes
        # find None values, and replace with 0
        for pk, priority in priorites.items():
            if priority is None:
                priorites[pk] = 0
                node = load_node(pk)
                node.base.extras.set("_scheduler_priority", 0)
                self.logger.report(f"Process {pk} has no priority, setting to 0")
        return priorites

    def reset(self) -> None:
        """Reset the scheduler node."""
        self.is_running = False
        self.waiting_process = []
        self.running_process = []
        self.running_workflow = []
        self.running_calcjob = []
        self.finished_process = []
        self.max_processes = 10000
        self.max_calcjobs = 100
        self.max_workflows = 100
        self.next_priority = 0
