from __future__ import annotations
from typing import Callable, Dict, Optional, Union
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida_workgraph.task import Task
from .workgraph import WorkGraph
import inspect
from .task import TaskHandle
from node_graph.node_spec import NodeSpec
from node_graph.socket_spec import SocketSpec
from aiida_workgraph.tasks.aiida import _build_aiida_function_nodespec
from node_graph.error_handler import ErrorHandlerSpec, normalize_error_handlers


def _spec_for(
    obj,
    *,
    identifier: Optional[str],
    inputs: Optional[SocketSpec] = None,
    outputs: Optional[SocketSpec] = None,
) -> NodeSpec:
    from aiida_workgraph.socket_spec import from_aiida_process
    from aiida_workgraph.utils import inspect_aiida_component_type
    from node_graph.executor import RuntimeExecutor

    # AiiDA process classes
    if inspect.isclass(obj) and issubclass(obj, (CalcJob, WorkChain)):
        from aiida_workgraph.tasks.aiida import AiiDAProcessTask

        in_spec, out_spec = from_aiida_process(obj)
        return NodeSpec(
            identifier=identifier or obj.__name__,
            mode='decorator_build',
            catalog='AIIDA',
            inputs=in_spec,
            outputs=out_spec,
            executor=RuntimeExecutor.from_callable(obj),
            base_class=AiiDAProcessTask,
            metadata={'node_type': inspect_aiida_component_type(obj)},
            decorator_path='aiida_workgraph.decorator.task',
        )

    # AiiDA process functions (calcfunction/workfunction)
    if callable(obj) and getattr(obj, 'node_class', False):
        from aiida_workgraph.tasks.aiida import _build_aiida_function_nodespec

        return _build_aiida_function_nodespec(obj, identifier=identifier, in_spec=inputs, out_spec=outputs)

    # Plain Python function -> PyFunction
    if callable(obj):
        from aiida_workgraph.tasks.pythonjob_tasks import build_pyfunction_nodespec

        return build_pyfunction_nodespec(obj, identifier=identifier, in_spec=inputs, out_spec=outputs)

    raise ValueError(f'Unsupported object for @task: {obj!r}')


def build_task_from_callable(
    executor: Callable,
    inputs: Optional[SocketSpec | list] = None,
    outputs: Optional[SocketSpec | list] = None,
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
        hasattr(executor, '_TaskCls')
        and inspect.isclass(executor._TaskCls)
        and issubclass(executor._TaskCls, Task)
        or inspect.isclass(executor)
        and issubclass(executor, Task)
    ):
        return executor
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, 'node_class', False):
            return task(inputs=inputs, outputs=outputs)(executor)
        else:
            return task(inputs=inputs, outputs=outputs)(executor)
    else:
        if issubclass(executor, CalcJob) or issubclass(executor, WorkChain):
            if inputs is not None or outputs is not None:
                raise ValueError('Can not override inputs or outputs of an AiiDA process classes.')
            return task()(executor)
    raise ValueError(f'The executor {executor} is not supported.')


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
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
        catalog: str = 'Others',
    ) -> Callable:
        """Generate a decorator that register a function as a task.

        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(obj: Union[WorkGraph, type, callable]) -> TaskHandle:
            spec = _spec_for(obj, identifier=identifier, inputs=inputs, outputs=outputs)
            handlers = normalize_error_handlers(error_handlers)
            # allow catalog override
            spec = NodeSpec(
                identifier=spec.identifier,
                mode=spec.mode,
                catalog=catalog or spec.catalog,
                inputs=spec.inputs,
                outputs=spec.outputs,
                executor=spec.executor,
                error_handlers=handlers,
                metadata=spec.metadata,
                base_class=spec.base_class,
                base_class_path=spec.base_class_path,
                decorator=spec.decorator,
                decorator_path=spec.decorator_path,
                version=spec.version,
            )
            handle = TaskHandle(spec)
            handle._func = obj
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def decorator_graph(
        identifier: Optional[str] = None,
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        max_depth: int = 100,
    ) -> Callable:
        """Generate a decorator that register a function as a graph task.
        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.graph_task import _build_graph_task_nodespec

            handle = TaskHandle(
                _build_graph_task_nodespec(
                    func,
                    identifier=identifier,
                    in_spec=inputs,
                    out_spec=outputs,
                    max_depth=max_depth,
                )
            )
            handle._func = func
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def calcfunction(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            func_decorated = calcfunction(func)
            handle = TaskHandle(
                _build_aiida_function_nodespec(
                    func_decorated,
                    in_spec=inputs,
                    out_spec=outputs,
                    error_handlers=error_handlers,
                )
            )
            handle._func = func_decorated
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def workfunction(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            func_decorated = workfunction(func)
            handle = TaskHandle(
                _build_aiida_function_nodespec(
                    func_decorated,
                    in_spec=inputs,
                    out_spec=outputs,
                    error_handlers=error_handlers,
                )
            )
            handle._func = func_decorated
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def pythonjob(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.pythonjob_tasks import build_pythonjob_nodespec

            spec = build_pythonjob_nodespec(
                func,
                in_spec=inputs,
                out_spec=outputs,
                error_handlers=error_handlers,
            )
            handle = TaskHandle(spec)
            handle._func = func
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def monitor(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.pythonjob_tasks import build_monitor_function_nodespec

            handle = TaskHandle(
                build_monitor_function_nodespec(
                    func,
                    in_spec=inputs,
                    out_spec=outputs,
                    error_handlers=error_handlers,
                )
            )
            handle._func = func
            return handle

        return decorator

    # Making decorator_task accessible as 'task'
    task = decorator_task

    # Making decorator_graph accessible as 'graph'
    graph = decorator_graph

    def __call__(self, *args, **kwargs):
        # This allows using '@task' to directly apply the decorator_task functionality
        if len(args) == 1 and isinstance(args[0], Callable) and len(kwargs) == 0:
            return self.decorator_task()(args[0])
        else:
            return self.decorator_task(*args, **kwargs)


task = TaskDecoratorCollection()
