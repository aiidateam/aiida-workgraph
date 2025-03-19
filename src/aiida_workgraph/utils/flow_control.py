from aiida_workgraph.socket import TaskSocket
from aiida_workgraph.tasks.task_pool import TaskPool

DEFAULT_MAP_PLACEHOLDER = "map_input"


class BaseFlowBlock:
    def __init__(self, condition: TaskSocket):
        self.condition = condition
        # The "work graph" is always the parent of the node in question
        self.wg = condition._parent._node.graph

    def _generate_name(self, prefix: str) -> str:
        # Subclasses can define their own identifier (e.g. "workgraph.if_zone" or "workgraph.while_zone")
        # for now, just do a dummy:
        existing = [
            task
            for task in self.wg.tasks
            if hasattr(task, "identifier") and task.identifier == self._task_identifier
        ]
        index = len(existing) + 1
        return f"{prefix}_{index}"

    @property
    def _task_identifier(self):
        """
        Subclasses must override to specify the identifier
        ("workgraph.if_zone" or "workgraph.while_zone", etc.)
        """
        raise NotImplementedError


class If_(BaseFlowBlock):
    _task_identifier = "workgraph.if_zone"

    def __init__(self, condition: TaskSocket):
        super().__init__(condition)
        self.true_zone = self.wg.add_task(
            TaskPool.workgraph.if_zone,
            name=self._generate_name("if_true"),
            conditions=self.condition,
        )
        self.false_zone = None

    def __call__(self, *tasks):
        _add_tasks_to_zone(self.true_zone, tasks)
        return self

    def else_(self, *tasks):
        self.false_zone = self.wg.add_task(
            TaskPool.workgraph.if_zone,
            name=self._generate_name("if_false"),
            conditions=self.condition,
            invert_condition=True,
        )
        _add_tasks_to_zone(self.false_zone, tasks)
        return self


class While_(BaseFlowBlock):
    _task_identifier = "workgraph.while_zone"

    def __init__(self, condition: TaskSocket, max_iterations: int = 10000):
        super().__init__(condition)
        self.zone = self.wg.add_task(
            TaskPool.workgraph.while_zone,
            name=self._generate_name("while"),
            conditions=self.condition,
            max_iterations=max_iterations,
        )

    def __call__(self, *tasks):
        _add_tasks_to_zone(self.zone, tasks)
        return self


class Map_(BaseFlowBlock):
    _task_identifier = "workgraph.map_zone"

    def __init__(self, source: TaskSocket, placeholder: str = DEFAULT_MAP_PLACEHOLDER):
        super().__init__(source)
        self.zone = self.wg.add_task(
            TaskPool.workgraph.map_zone,
            name=self._generate_name("map"),
            source=self.condition,
            placeholder=placeholder,
        )

    def __call__(self, *tasks):
        _add_tasks_to_zone(self.zone, tasks)
        return self


def _add_tasks_to_zone(zone, tasks):
    """
    A helper that takes a 'zone' (i.e. the node to which tasks will be attached)
    and a list of tasks (some of which might be normal tasks, while_ objects, or
    if_ objects). It dispatches them to the zone accordingly.
    """
    from aiida_workgraph.task import Task
    from node_graph.socket import NodeSocket

    normal_tasks = []
    for t in tasks:
        if isinstance(t, Task):
            normal_tasks.append(t)
        elif isinstance(t, NodeSocket):
            normal_tasks.append(t._node)
    zone.children.add(normal_tasks)

    while_tasks = [t for t in tasks if isinstance(t, While_)]
    for w_task in while_tasks:
        zone.children.add(w_task.zone)

    if_tasks = [t for t in tasks if isinstance(t, If_)]
    for i_task in if_tasks:
        zone.children.add(i_task.true_zone)
        if i_task.false_zone:
            zone.children.add(i_task.false_zone)


def if_(condition):
    return If_(condition)


def while_(condition, max_iterations=10000):
    return While_(condition, max_iterations)


def map_(source, placeholder: str = DEFAULT_MAP_PLACEHOLDER):
    return Map_(source, placeholder)


map_.default_placeholder = DEFAULT_MAP_PLACEHOLDER
