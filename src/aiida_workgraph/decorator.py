from __future__ import annotations
import functools
from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida_workgraph.task import Task
from .workgraph import WorkGraph
import inspect
from node_graph.executor import NodeExecutor
from aiida_workgraph.tasks.factory import (
    DecoratedFunctionTaskFactory,
    AiiDAComponentTaskFactory,
    WorkGraphTaskFactory,
    PyFunctionTaskFactory,
)


def build_task(
    executor: Union[Callable, str],
    inputs: Optional[List[str | dict]] = None,
    outputs: Optional[List[str | dict]] = None,
) -> Task:
    """Build task from executor."""

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
        hasattr(executor, "_TaskCls")
        and inspect.isclass(executor._TaskCls)
        and issubclass(executor._TaskCls, Task)
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
            return PyFunctionTaskFactory.from_function(
                executor,
                inputs=inputs,
                outputs=outputs,
            )
    else:
        if issubclass(executor, CalcJob) or issubclass(executor, WorkChain):
            return AiiDAComponentTaskFactory.from_aiida_component(
                executor, inputs=inputs, outputs=outputs
            )
    raise ValueError(f"The executor {executor} is not supported.")


def _assign_wg_outputs(
    outputs: Any, wg: WorkGraph, graph_task_output_names: List[str]
) -> None:
    """
    Inspect the raw outputs from the function and attach them to the WorkGraph.
    """
    from node_graph.socket import BaseSocket

    if isinstance(outputs, BaseSocket):
        wg.outputs.result = outputs
    elif isinstance(outputs, dict):
        wg.outputs = outputs
    elif isinstance(outputs, tuple):
        if len(outputs) != len(graph_task_output_names):
            raise ValueError(
                f"The length of the outputs {len(outputs)} does not match the length of the \
                    Graph task outputs {len(graph_task_output_names)}."
            )
        outputs_dict = {}
        for i, output in enumerate(outputs):
            outputs_dict[graph_task_output_names[i]] = output
        wg.outputs = outputs_dict
    else:
        wg.outputs.result = outputs


def _run_func_with_wg(
    func: Callable,
    graph_task_output_names: List[str],
    args: tuple,
    kwargs: dict,
    var_kwargs: Optional[dict] = None,
) -> WorkGraph:
    """
    Run func(*args, **kwargs, **(var_kwargs or {})) inside a WorkGraph,
    assign its outputs and return the WorkGraph.
    """
    merged = {**kwargs, **(var_kwargs or {})}
    with WorkGraph(func.__name__) as wg:
        raw = func(*args, **merged)
        _assign_wg_outputs(raw, wg, graph_task_output_names)
        return wg


def _make_wrapper(TaskCls, func=None):
    """
    Common wrapper that, when called, adds a node to the current graph
    and returns the outputs.
    """

    @functools.wraps(func)
    def wrapper(*call_args, **call_kwargs):
        from aiida_workgraph.manager import get_current_graph

        graph = get_current_graph()
        if graph is None:
            raise RuntimeError(f"No active Graph available for {TaskCls.name}.")
        task = graph.add_task(TaskCls)
        active_zone = getattr(graph, "_active_zone", None)
        if active_zone:
            active_zone.children.add(task)

        inputs = dict(call_kwargs or {})
        if func is not None:
            arguments = list(call_args)
            orginal_func = func._func if hasattr(func, "_func") else func

            for name, parameter in inspect.signature(orginal_func).parameters.items():
                if parameter.kind in [
                    parameter.POSITIONAL_ONLY,
                    parameter.POSITIONAL_OR_KEYWORD,
                ]:
                    try:
                        inputs[name] = arguments.pop(0)
                    except IndexError:
                        pass
                elif parameter.kind is parameter.VAR_POSITIONAL:
                    # not supported
                    raise ValueError("VAR_POSITIONAL is not supported.")
        task.set(inputs)
        return task.outputs

    # Expose the TaskCls on the wrapper if you want
    wrapper._TaskCls = wrapper._NodeCls = TaskCls
    wrapper._func = func
    wrapper.is_decoratored = True
    return wrapper


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
        store_provenance: bool = True,
    ) -> Callable:
        """Generate a decorator that register a function as a task.

        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            properties (list): task properties
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(callable):

            # function or builtin function
            if inspect.isfunction(callable) or callable.__module__ == "builtins":
                # calcfunction and workfunction
                if getattr(callable, "node_class", False):
                    TaskCls = AiiDAComponentTaskFactory.from_aiida_component(
                        callable, inputs=inputs, outputs=outputs
                    )
                else:
                    if store_provenance:
                        TaskCls = PyFunctionTaskFactory.from_function(
                            callable,
                            inputs=inputs,
                            outputs=outputs,
                            error_handlers=error_handlers,
                        )
                    else:
                        TaskCls = DecoratedFunctionTaskFactory.from_function(
                            func=callable,
                            identifier=identifier,
                            task_type=task_type,
                            properties=properties,
                            inputs=inputs,
                            outputs=outputs,
                            error_handlers=error_handlers,
                            catalog=catalog,
                        )
            else:
                if issubclass(callable, CalcJob) or issubclass(callable, WorkChain):
                    TaskCls = AiiDAComponentTaskFactory.from_aiida_component(
                        callable, inputs=inputs, outputs=outputs
                    )
            # if callable is a function, we pass it to the make_wrapper
            if not inspect.isfunction(callable):
                callable = None
            return _make_wrapper(TaskCls, func=callable)

        return decorator

    @staticmethod
    @nonfunctional_usage
    def decorator_graph(
        identifier: Optional[str] = None,
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Tuple[str, str]]] = None,
        outputs: Optional[List[Tuple[str, str]]] = None,
        catalog: str = "Others",
    ) -> Callable:
        """Generate a decorator that register a function as a graph task.
        Attributes:
            indentifier (str): task identifier
            catalog (str): task catalog
            properties (list): task properties
            inputs (list): task inputs
            outputs (list): task outputs
        """

        def decorator(func):
            from aiida_workgraph.tasks.builtins import GraphTask

            TaskCls = DecoratedFunctionTaskFactory.from_function(
                func=func,
                identifier=identifier,
                task_type="graph_task",
                properties=properties,
                inputs=inputs,
                outputs=outputs,
                catalog=catalog,
                node_class=GraphTask,
            )

            wrapped_func = _make_wrapper(TaskCls, func)

            def build_graph(*args, **kwargs):
                """This function is used to get the graph from the wrapped function."""
                graph_task_output_names = [
                    name
                    for name, socket in TaskCls._ndata.get("outputs", {})
                    .get("sockets", {})
                    .items()
                    if not socket.get("metadata", {}).get("builtin_socket", False)
                ]

                return _run_func_with_wg(func, graph_task_output_names, args, kwargs)

            wrapped_func.build_graph = build_graph

            return wrapped_func

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

            return _make_wrapper(TaskCls, func_decorated)

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

            return _make_wrapper(TaskCls, func_decorated)

        return decorator

    @staticmethod
    @nonfunctional_usage
    def pythonjob(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
    ) -> Callable:
        def decorator(func):
            from aiida_workgraph.tasks.factory import (
                PythonJobTaskFactory,
            )

            TaskCls = PythonJobTaskFactory.from_function(
                func, inputs=inputs, outputs=outputs, error_handlers=error_handlers
            )

            return _make_wrapper(TaskCls, func)

        return decorator

    @staticmethod
    @nonfunctional_usage
    def awaitable(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
    ) -> Callable:
        def decorator(func):
            from aiida_workgraph.tasks.factory.awaitable_task import (
                AwaitableFunctionTaskFactory,
            )

            TaskCls = AwaitableFunctionTaskFactory.from_function(
                func=func,
                inputs=inputs,
                outputs=outputs,
            )

            return _make_wrapper(TaskCls, func)

        return decorator

    @staticmethod
    @nonfunctional_usage
    def monitor(
        inputs: Optional[List[str | dict]] = None,
        outputs: Optional[List[str | dict]] = None,
    ) -> Callable:
        def decorator(func):
            from aiida_workgraph.tasks.factory.awaitable_task import (
                MonitorFunctionTaskFactory,
            )

            sig = inspect.signature(func)
            forbidden = {"interval", "timeout"}
            conflicts = forbidden.intersection(sig.parameters)
            if conflicts:
                raise ValueError(
                    f"Function '{func.__name__}' defines parameter(s) {sorted(conflicts)}, "
                    "which conflict with default 'monitor' arguments."
                )
            task_inputs = inputs if inputs is not None else {}
            # Add default interval and timeout
            task_inputs.update(
                {
                    "interval": {
                        "identifier": "workgraph.float",
                        "property": {"default": 5},
                    },
                    "timeout": {
                        "identifier": "workgraph.float",
                        "property": {"default": 3600},
                    },
                }
            )
            TaskCls = MonitorFunctionTaskFactory.from_function(
                func=func,
                inputs=task_inputs,
                outputs=outputs,
            )

            return _make_wrapper(TaskCls, func)

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
