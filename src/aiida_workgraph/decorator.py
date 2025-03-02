from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida_workgraph.task import Task
import inspect
from node_graph.executor import NodeExecutor
from aiida_workgraph.tasks.factory import (
    DecoratedFunctionTaskFactory,
    AiiDAComponentTaskFactory,
    WorkGraphTaskFactory,
)


def build_task(
    executor: Union[Callable, str],
    inputs: Optional[List[str | dict]] = None,
    outputs: Optional[List[str | dict]] = None,
) -> Task:
    """Build task from executor."""
    from aiida_workgraph.workgraph import WorkGraph

    if isinstance(executor, WorkGraph):
        return WorkGraphTaskFactory.create_task(executor)
    elif isinstance(executor, str):
        executor = NodeExecutor(module_path=executor).executor
    if callable(executor):
        return build_task_from_callable(executor, inputs=inputs, outputs=outputs)


def build_task_from_callable(
    executor: Callable,
    inputs: Optional[List[str | dict]] = None,
    outputs: Optional[List[str | dict]] = None,
) -> Task:
    """Build task from a callable object.
    First, check if the executor is already a task.
    If not, check if it is a function or a class.
    If it is a function, build task from function.
    If it is a class, it only supports CalcJob and WorkChain.
    """
    from aiida_workgraph.task import Task

    # if it is already a task, return it
    if (
        hasattr(executor, "TaskCls")
        and inspect.isclass(executor.TaskCls)
        and issubclass(executor.TaskCls, Task)
        or inspect.isclass(executor)
        and issubclass(executor, Task)
    ):
        return executor
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, "node_class", False):
            return AiiDAComponentTaskFactory.from_aiida_component(
                executor, inputs=inputs, outputs=outputs
            )
        else:
            return DecoratedFunctionTaskFactory.from_function(
                executor, inputs=inputs, outputs=outputs
            )
    else:
        if issubclass(executor, CalcJob):
            return AiiDAComponentTaskFactory.from_aiida_component(
                executor, inputs=inputs, outputs=outputs
            )
        elif issubclass(executor, WorkChain):
            return AiiDAComponentTaskFactory.from_aiida_component(
                executor, inputs=inputs, outputs=outputs
            )
    raise ValueError(f"The executor {executor} is not supported.")


def nonfunctional_usage(callable: Callable):
    """
    This is a decorator for a decorator factory (a function that returns a decorator).
    It allows the usage of the decorator factory in a nonfunctional way. So a decorator
    factory that has been decorated by this decorator that could only be used befor like
    this

    .. code-block:: python

        @decorator_factory()
        def foo():
            pass

    can now be also used like this

    .. code-block:: python

        @decorator_factory
        def foo():
            pass

    """

    def decorator_task_wrapper(*args, **kwargs):
        if len(args) == 1 and isinstance(args[0], Callable) and len(kwargs) == 0:
            return callable()(args[0])
        else:
            return callable(*args, **kwargs)

    return decorator_task_wrapper


class TaskDecoratorCollection:
    """Collection of task decorators."""

    @staticmethod
    @nonfunctional_usage
    def decorator_task(
        identifier: Optional[str] = None,
        task_type: str = "Normal",
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
    ) -> Callable:
        """Generate a decorator that register a function as a task.

        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            properties (list): task properties
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(func):
            # at least one output is required
            task_outputs = outputs or [
                {"identifier": "workgraph.any", "name": "result"}
            ]
            TaskCls = DecoratedFunctionTaskFactory.from_function(
                func=func,
                identifier=identifier,
                task_type=task_type,
                properties=properties,
                inputs=inputs,
                outputs=task_outputs,
                error_handlers=error_handlers,
                catalog=catalog,
            )

            func.identifier = TaskCls.identifier
            func.TaskCls = func.NodeCls = TaskCls
            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def decorator_graph_builder(
        identifier: Optional[str] = None,
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Tuple[str, str]]] = None,
        outputs: Optional[List[Tuple[str, str]]] = None,
        catalog: str = "Others",
    ) -> Callable:
        """Generate a decorator that register a function as a graph_builder task.
        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            properties (list): task properties
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(func):
            from aiida_workgraph.tasks.builtins import GraphBuilderTask

            task_outputs = [
                {"identifier": "workgraph.any", "name": output["name"]}
                for output in outputs or []
            ]
            TaskCls = DecoratedFunctionTaskFactory.from_function(
                func=func,
                identifier=identifier,
                task_type="graph_builder",
                properties=properties,
                inputs=inputs,
                outputs=task_outputs,
                catalog=catalog,
                group_inputs=inputs,
                group_outputs=outputs,
                node_class=GraphBuilderTask,
            )
            func.TaskCls = func.NodeCls = TaskCls
            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def calcfunction(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
    ) -> Callable:
        def decorator(func):
            func_decorated = calcfunction(func)
            TaskCls = AiiDAComponentTaskFactory.from_aiida_component(
                func_decorated,
                inputs=inputs,
                outputs=outputs,
                error_handlers=error_handlers,
            )
            func_decorated.identifier = TaskCls.identifier
            func_decorated.TaskCls = func_decorated.NodeCls = TaskCls
            return func_decorated

        return decorator

    @staticmethod
    @nonfunctional_usage
    def workfunction(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
    ) -> Callable:
        def decorator(func):
            func_decorated = workfunction(func)
            TaskCls = AiiDAComponentTaskFactory.from_aiida_component(
                func_decorated,
                inputs=inputs,
                outputs=outputs,
                error_handlers=error_handlers,
            )
            func_decorated.identifier = TaskCls.identifier
            func_decorated.TaskCls = func_decorated.NodeCls = TaskCls

            return func_decorated

        return decorator

    @staticmethod
    @nonfunctional_usage
    def pythonjob(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
    ) -> Callable:
        def decorator(func):
            from aiida_workgraph.tasks.factory.pythonjob_task import (
                PythonJobTaskFactory,
            )

            TaskCls = PythonJobTaskFactory.from_function(
                func, inputs=inputs, outputs=outputs, error_handlers=error_handlers
            )
            func.TaskCls = func.NodeCls = TaskCls

            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def awaitable(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
    ) -> Callable:
        def decorator(func):
            # at least one output is required
            task_outputs = outputs or [
                {"identifier": "workgraph.any", "name": "result"}
            ]
            TaskCls = DecoratedFunctionTaskFactory.from_function(
                func=func,
                task_type="awaitable",
                inputs=inputs,
                outputs=task_outputs,
            )
            func.TaskCls = func.NodeCls = TaskCls
            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def monitor(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
    ) -> Callable:
        def decorator(func):
            # at least one output is required
            task_outputs = outputs or [
                {"identifier": "workgraph.any", "name": "result"}
            ]
            TaskCls = DecoratedFunctionTaskFactory.from_function(
                func=func,
                task_type="monitor",
                inputs=inputs,
                outputs=task_outputs,
            )
            # add an input: interval
            func.TaskCls = func.NodeCls = TaskCls
            return func

        return decorator

    # Making decorator_task accessible as 'task'
    task = decorator_task

    # Making decorator_graph_builder accessible as 'graph_builder'
    graph_builder = decorator_graph_builder

    def __call__(self, *args, **kwargs):
        # This allows using '@task' to directly apply the decorator_task functionality
        if len(args) == 1 and isinstance(args[0], Callable) and len(kwargs) == 0:
            return self.decorator_task()(args[0])
        else:
            return self.decorator_task(*args, **kwargs)


task = TaskDecoratorCollection()
