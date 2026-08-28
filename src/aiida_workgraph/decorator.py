from __future__ import annotations
from dataclasses import replace
from typing import Callable, Dict, Optional, Type, Union
from pydantic import BaseModel
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida_workgraph.task import Task
from .workgraph import WorkGraph
import inspect
from .task import TaskHandle
from node_graph.task_spec import TaskSpec
from node_graph.socket_spec import SocketSpec
from aiida_workgraph.socket_spec import SocketSpecAPI, node_typed_paths
from aiida_workgraph.tasks.aiida import _build_aiida_function_taskspec
from node_graph.error_handler import ErrorHandlerSpec, normalize_error_handlers
from aiida_workgraph.tasks.pythonjob_tasks import build_pyfunction_taskspec
from aiida_workgraph.tasks.aiida import AiiDAProcessTask
from node_graph.executor import RuntimeExecutor
from node_graph.input_model import ModelContractError, apply_models, rebind_executor_callable


def _spec_for(
    obj,
    *,
    identifier: Optional[str],
    inputs: Optional[SocketSpec] = None,
    outputs: Optional[SocketSpec] = None,
    catalog: str = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> TaskSpec:
    # AiiDA process classes
    if inspect.isclass(obj) and issubclass(obj, (CalcJob, WorkChain)):
        return AiiDAProcessTask.build(obj, attached_error_handlers=error_handlers)

    # AiiDA process functions (calcfunction/workfunction)
    if callable(obj) and getattr(obj, 'node_class', False):
        return _build_aiida_function_taskspec(
            obj,
            identifier=identifier,
            in_spec=inputs,
            out_spec=outputs,
            error_handlers=error_handlers,
            catalog=catalog or 'Others',
        )

    # Plain Python function -> PyFunction
    if callable(obj):
        spec = build_pyfunction_taskspec(
            obj,
            identifier=identifier,
            in_spec=inputs,
            out_spec=outputs,
            error_handlers=error_handlers,
            catalog=catalog or 'Others',
        )
        return spec

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
    from node_graph.task import Task

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


def _refuse_node_typed_inputs(model: Optional[Type[BaseModel]], spec: Optional[SocketSpec]) -> None:
    """Raise when a model asks a PyFunction body for a value it cannot be handed.

    A PyFunction's inputs are read out of their nodes before its body runs, so
    a field declaring an AiiDA type would be handed what the node carries and
    the model would refuse the very thing it asked for. A calcfunction's body
    is handed the nodes themselves, and is where such a field belongs.
    """
    if model is None or spec is None:
        return
    paths = node_typed_paths(spec)
    if not paths:
        return
    listed = ', '.join(repr(path) for path in paths)
    raise ModelContractError(
        f'{model.__name__} declares {listed} as an AiiDA type, and a task declared with '
        '@task runs its body as a PyFunction, which is handed the value a node carries, '
        'never the node.\n'
        'How to fix: declare the task with @task.calcfunction, whose body is handed the '
        'node; or declare the field as the Python type the body reads.'
    )


def _refuse_namespace_inputs(model: Optional[Type[BaseModel]], spec: Optional[SocketSpec]) -> None:
    """Raise when a model gives a calcfunction a parameter AiiDA cannot express.

    A process function's parameter is one port carrying one node. A field
    declaring a nested model, or a ``dict[str, T]``, asks for a namespace, and
    AiiDA turns the mapping it is handed into a single ``orm.Dict``: the
    members lose the nodes they were, and the model refuses the ``Dict`` it
    never declared.
    """
    if model is None or spec is None:
        return
    namespaces = [name for name, field in (spec.fields or {}).items() if field.is_namespace()]
    if not namespaces:
        return
    listed = ', '.join(repr(name) for name in namespaces)
    raise ModelContractError(
        f'{model.__name__} declares {listed} as a namespace -- a nested model or a '
        'dict[str, T] -- and a calcfunction parameter is one port carrying one node.\n'
        'How to fix: declare the task with @task, whose body is handed a namespace as a '
        'mapping; or give the model one field per value the body reads.'
    )


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
        input_model: Optional[Type[BaseModel]] = None,
        output_model: Optional[Type[BaseModel]] = None,
    ) -> Callable:
        """Generate a decorator that register a function as a task.

        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            inputs (list): task inputs
            outputs (list): task outputs
            input_model (BaseModel): model declaring the input sockets, checked at every
                call and again before the body runs
            output_model (BaseModel): model declaring the output sockets and validating
                the return value
        """

        def decorator(obj: Union[WorkGraph, type, callable]) -> TaskHandle:
            normalized_handlers = normalize_error_handlers(error_handlers)
            in_spec, out_spec, executor = apply_models(
                obj, inputs, outputs, input_model, output_model, api=SocketSpecAPI
            )
            _refuse_node_typed_inputs(input_model, in_spec)
            spec = _spec_for(
                obj,
                identifier=identifier,
                catalog=catalog,
                inputs=in_spec,
                outputs=out_spec,
                error_handlers=normalized_handlers,
            )
            if executor is not obj:
                # The spec is inferred from the undecorated function, so its
                # signature, source and return annotation stay visible; only
                # what runs changes.
                spec = replace(spec, executor=RuntimeExecutor.from_callable(executor))

            handle = TaskHandle(spec)
            handle._callable = executor
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def decorator_graph(
        identifier: Optional[str] = None,
        catalog: Optional[str] = None,
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        max_depth: int = 100,
        max_number_jobs: Optional[int] = None,
        input_model: Optional[Type[BaseModel]] = None,
        output_model: Optional[Type[BaseModel]] = None,
    ) -> Callable:
        """Generate a decorator that register a function as a graph task.
        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            inputs (list): task inputs
            outputs (list): task outputs
            input_model (BaseModel): model declaring the input sockets, checked at every
                call and again when the graph is expanded
            output_model (BaseModel): refused; a graph returns socket references, which
                stand for values that do not exist yet
        """

        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.graph_task import _build_graph_task_taskspec

            in_spec, _, executor = apply_models(
                func, inputs, None, input_model, output_model, is_graph=True, api=SocketSpecAPI
            )
            handle = TaskHandle(
                _build_graph_task_taskspec(
                    func,
                    identifier=identifier,
                    catalog=catalog,
                    in_spec=in_spec,
                    out_spec=outputs,
                    max_depth=max_depth,
                    max_number_jobs=max_number_jobs,
                )
            )
            handle._callable = executor
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def calcfunction(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        catalog: Optional[str] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
        input_model: Optional[Type[BaseModel]] = None,
        output_model: Optional[Type[BaseModel]] = None,
    ) -> Callable:
        """Generate a decorator registering a function as a calcfunction task.

        Attributes:
            inputs (list): task inputs
            outputs (list): task outputs
            input_model (BaseModel): model declaring the input sockets, checked at every
                call and again before the body runs; a field declaring an AiiDA type is
                handed the node, which is what a calcfunction's body receives
            output_model (BaseModel): model declaring the output sockets and validating
                the return value
        """

        def decorator(func) -> TaskHandle:
            in_spec, out_spec, executor = apply_models(
                func, inputs, outputs, input_model, output_model, api=SocketSpecAPI
            )
            _refuse_namespace_inputs(input_model, in_spec)
            # The models are enforced inside the process, so what AiiDA runs is
            # the wrapper and the nodes it is called with reach the body. The
            # calcfunction is what the executor has to resolve to, so it takes
            # over the name the wrapper was bound under.
            func_decorated = calcfunction(executor)
            rebind_executor_callable(func_decorated, executor)
            handle = TaskHandle(
                _build_aiida_function_taskspec(
                    func_decorated,
                    in_spec=in_spec,
                    out_spec=out_spec,
                    catalog=catalog,
                    error_handlers=error_handlers,
                )
            )
            handle._callable = func_decorated
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def workfunction(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        catalog: Optional[str] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            func_decorated = workfunction(func)
            handle = TaskHandle(
                _build_aiida_function_taskspec(
                    func_decorated,
                    in_spec=inputs,
                    out_spec=outputs,
                    catalog=catalog,
                    error_handlers=error_handlers,
                )
            )
            handle._callable = func_decorated
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def pythonjob(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        catalog: Optional[str] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.pythonjob_tasks import build_pythonjob_taskspec

            spec = build_pythonjob_taskspec(
                func,
                in_spec=inputs,
                out_spec=outputs,
                catalog=catalog,
                error_handlers=error_handlers,
            )
            handle = TaskHandle(spec)
            handle._callable = func
            return handle

        return decorator

    @staticmethod
    @nonfunctional_usage
    def monitor(
        inputs: Optional[SocketSpec | list] = None,
        outputs: Optional[SocketSpec | list] = None,
        catalog: Optional[str] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ) -> Callable:
        def decorator(func) -> TaskHandle:
            from aiida_workgraph.tasks.pythonjob_tasks import build_monitor_function_taskspec

            handle = TaskHandle(
                build_monitor_function_taskspec(
                    func,
                    in_spec=inputs,
                    out_spec=outputs,
                    catalog=catalog,
                    error_handlers=error_handlers,
                )
            )
            handle._callable = func
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
