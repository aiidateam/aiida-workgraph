from aiida_workgraph.socket import TaskSocket
from aiida_workgraph.tasks.task_pool import TaskPool


def if_(condition):
    """Helper function to create an If_ object."""
    return If_(condition)


class If_:
    def __init__(self, condition: TaskSocket):
        self.condition = condition
        self.wg = condition._parent._node.parent
        self.true_zone = self.wg.add_task(
            TaskPool.workgraph.if_zone, name=self._generate_name(), conditions=condition
        )
        self.false_zone = None

    def _generate_name(self, prefix: str = "if_true") -> str:
        n = (
            len(
                [
                    task
                    for task in self.wg.tasks
                    if task.identifier == "workgraph.if_zone"
                ]
            )
            + 1
        )
        return f"{prefix}_{n}"

    def __call__(self, *tasks):
        """
        Called when the user does:
            if_(condition)(<tasks-for-true>)
        """
        from aiida_workgraph.task import Task

        tasks = [task for task in tasks if isinstance(task, Task)]
        self.true_zone.children.add(tasks)
        return self

    def elif_(self):
        """Not implemented."""
        raise NotImplementedError("elif_ not implemented.")

    def else_(self, *tasks):
        """
        Called when the user does:
            .else_(<tasks-for-false>)
        """
        self.false_zone = self.wg.add_task(
            TaskPool.workgraph.if_zone,
            name=self._generate_name("if_false"),
            conditions=self.condition,
            invert_condition=True,
        )
        from aiida_workgraph.task import Task

        tasks = [task for task in tasks if isinstance(task, Task)]
        self.false_zone.children.add(tasks)

        return self


def while_(condition, max_iterations: int = 10000):
    """Helper function to create an While_ object."""
    return While_(condition, max_iterations=max_iterations)


class While_:
    def __init__(self, condition: TaskSocket, max_iterations: int = 10000):
        self.condition = condition
        self.wg = condition._parent._node.parent
        self.zone = self.wg.add_task(
            TaskPool.workgraph.while_zone,
            name=self._generate_name(),
            conditions=condition,
            max_iterations=max_iterations,
        )

    def _generate_name(self, prefix: str = "while") -> str:
        n = (
            len(
                [
                    task
                    for task in self.wg.tasks
                    if task.identifier == "workgraph.while_zone"
                ]
            )
            + 1
        )
        return f"{prefix}_{n}"

    def __call__(self, *tasks):
        """
        Called when the user does:
            while_(condition)(<tasks-for-loop>)
        """
        from aiida_workgraph.task import Task

        tasks = [task for task in tasks if isinstance(task, Task)]
        self.zone.children.add(tasks)
        return self
