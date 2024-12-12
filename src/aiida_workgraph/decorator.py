from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida_workgraph.utils import get_executor
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida_workgraph.task import Task
from aiida_workgraph.utils import build_callable, validate_task_inout
import inspect
from aiida_workgraph.config import builtin_inputs, builtin_outputs, task_types
from aiida_workgraph.orm.mapping import type_mapping


def create_task(tdata):
    """Wrap create_node from node_graph to create a Task."""
    from node_graph.decorator import create_node
    from node_graph.utils import list_to_dict

    tdata["type_mapping"] = type_mapping
    tdata["metadata"]["node_type"] = tdata["metadata"].pop("task_type")
    tdata["properties"] = list_to_dict(tdata.get("properties", {}))
    tdata["inputs"] = list_to_dict(tdata.get("inputs", {}))
    tdata["outputs"] = list_to_dict(tdata.get("outputs", {}))

    return create_node(tdata)


def build_task(
    executor: Union[Callable, str],
    inputs: Optional[List[str | dict]] = None,
    outputs: Optional[List[str | dict]] = None,
) -> Task:
    """Build task from executor."""
    from aiida_workgraph.workgraph import WorkGraph

    if inputs:
        inputs = validate_task_inout(inputs, "inputs")

    if outputs:
        outputs = validate_task_inout(outputs, "outputs")

    if isinstance(executor, WorkGraph):
        return build_task_from_workgraph(executor)
    elif isinstance(executor, str):
        (
            module,
            executor_name,
        ) = executor.rsplit(".", 1)
        executor, _ = get_executor(
            {"module_path": module, "callable_name": executor_name}
        )
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
        hasattr(executor, "task")
        and inspect.isclass(executor.task)
        and issubclass(executor.task, Task)
        or inspect.isclass(executor)
        and issubclass(executor, Task)
    ):
        return executor
    tdata = {"metadata": {}}
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, "node_class", False):
            tdata["metadata"]["task_type"] = task_types.get(
                executor.node_class, "NORMAL"
            )
            tdata["callable"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        else:
            tdata["metadata"]["task_type"] = "NORMAL"
            tdata["executor"] = executor
            return build_task_from_function(executor, inputs=inputs, outputs=outputs)
    else:
        if issubclass(executor, CalcJob):
            tdata["metadata"]["task_type"] = "CALCJOB"
            tdata["callable"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        elif issubclass(executor, WorkChain):
            tdata["metadata"]["task_type"] = "WORKCHAIN"
            tdata["callable"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
    raise ValueError(f"The executor {executor} is not supported.")


def build_task_from_function(
    executor: Callable,
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Task:
    """Build task from function."""
    return TaskDecoratorCollection.decorator_task(inputs=inputs, outputs=outputs)(
        executor
    ).task


def build_task_from_AiiDA(
    tdata: Dict[str, Any],
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Tuple[Task, Dict[str, Any]]:
    """Register a task from a AiiDA component.
    For example: CalcJob, WorkChain, CalcFunction, WorkFunction."""
    from aiida_workgraph.utils.inspect_aiida_components import (
        get_task_data_from_aiida_component,
    )

    tdata = get_task_data_from_aiida_component(
        tdata=tdata, inputs=inputs, outputs=outputs
    )
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata


def build_pythonjob_task(func: Callable) -> Task:
    """Build PythonJob task from function."""
    from aiida_pythonjob import PythonJob
    from copy import deepcopy

    # if the function is not a task, build a task from the function
    if not hasattr(func, "node"):
        TaskDecoratorCollection.decorator_task()(func)
    if func.node.node_type.upper() == "GRAPH_BUILDER":
        raise ValueError(
            "GraphBuilder task cannot be run remotely. Please remove 'PythonJob'."
        )
    tdata = {
        "metadata": {"task_type": "PYTHONJOB"},
        "callable": PythonJob,
    }
    _, tdata_py = build_task_from_AiiDA(tdata)
    tdata = deepcopy(func.tdata)
    # merge the inputs and outputs from the PythonJob task to the function task
    # skip the already existed inputs and outputs
    for input in [
        {"identifier": "workgraph.string", "name": "computer"},
        {"identifier": "workgraph.any", "name": "command_info"},
    ]:
        input["list_index"] = len(tdata["inputs"]) + 1
        tdata["inputs"][input["name"]] = input
    for name, input in tdata_py["inputs"].items():
        if name not in tdata["inputs"]:
            input["list_index"] = len(tdata["inputs"]) + 1
            tdata["inputs"][name] = input
    for name, output in tdata_py["outputs"].items():
        if name not in tdata["outputs"]:
            output["list_index"] = len(tdata["outputs"]) + 1
            tdata["outputs"][name] = output
    for output in [{"identifier": "workgraph.any", "name": "exit_code"}]:
        output["list_index"] = len(tdata["outputs"]) + 1
        tdata["outputs"][output["name"]] = output
    # change "copy_files" link_limit to 1e6
    tdata["inputs"]["copy_files"]["link_limit"] = 1e6
    tdata["metadata"]["task_type"] = "PYTHONJOB"
    tdata["identifier"] = "workgraph.pythonjob"
    tdata["metadata"]["node_class"] = {
        "module_path": "aiida_workgraph.tasks.pythonjob",
        "callable_name": "PythonJob",
    }
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata


def build_shelljob_task(outputs: list = None, parser_outputs: list = None) -> Task:
    """Build ShellJob with custom inputs and outputs."""
    from aiida_shell.calculations.shell import ShellJob
    from aiida_shell.parsers.shell import ShellParser

    tdata = {
        "metadata": {"task_type": "SHELLJOB"},
        "callable": ShellJob,
    }
    _, tdata = build_task_from_AiiDA(tdata)
    # Extend the outputs
    for output in [
        {"identifier": "workgraph.any", "name": "stdout"},
        {"identifier": "workgraph.any", "name": "stderr"},
    ]:
        output["list_index"] = len(tdata["outputs"]) + 1
        tdata["outputs"][output["name"]] = output
    outputs = [] if outputs is None else outputs
    parser_outputs = [] if parser_outputs is None else parser_outputs
    parser_outputs = validate_task_inout(parser_outputs, "parser_outputs")

    outputs = [
        {"identifier": "workgraph.any", "name": ShellParser.format_link_label(output)}
        for output in outputs
    ]
    outputs.extend(parser_outputs)
    # add user defined outputs
    for output in outputs:
        if output["name"] not in tdata["outputs"]:
            output["list_index"] = len(tdata["outputs"]) + 1
            tdata["outputs"][output["name"]] = output
    #
    tdata["identifier"] = "ShellJob"
    for input in [
        {"identifier": "workgraph.any", "name": "command"},
        {"identifier": "workgraph.any", "name": "resolve_command"},
    ]:
        input["list_index"] = len(tdata["inputs"]) + 1
        tdata["inputs"][input["name"]] = input
    tdata["metadata"]["task_type"] = "SHELLJOB"
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata


def build_task_from_workgraph(wg: any) -> Task:
    """Build task from workgraph."""

    tdata = {"metadata": {"task_type": "workgraph"}}
    inputs = []
    outputs = []
    group_outputs = []
    # add all the inputs/outputs from the tasks in the workgraph
    builtin_input_names = [input["name"] for input in builtin_inputs]
    builtin_output_names = [output["name"] for output in builtin_outputs]
    for task in wg.tasks:
        # inputs
        inputs.append(
            {
                "identifier": "workgraph.namespace",
                "name": f"{task.name}",
            }
        )
        for socket in task.inputs:
            if socket._name in builtin_input_names:
                continue
            inputs.append(
                {
                    "identifier": socket._identifier,
                    "name": f"{task.name}.{socket._name}",
                }
            )
        # outputs
        outputs.append(
            {
                "identifier": "workgraph.namespace",
                "name": f"{task.name}",
            }
        )
        for socket in task.outputs:
            if socket._name in builtin_output_names:
                continue
            outputs.append(
                {
                    "identifier": socket._identifier,
                    "name": f"{task.name}.{socket._name}",
                }
            )
            group_outputs.append(
                {
                    "name": f"{task.name}.{socket._name}",
                    "from": f"{task.name}.{socket._name}",
                }
            )
    # add built-in sockets
    for output in builtin_outputs:
        outputs.append(output.copy())
    for input in builtin_inputs:
        inputs.append(input.copy())
    tdata["metadata"]["node_class"] = {
        "module_path": "aiida_workgraph.task",
        "callable_name": "Task",
    }
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    tdata["identifier"] = wg.name
    # get wgdata from the workgraph
    wgdata = wg.prepare_inputs()["wg"]
    executor = {
        "module_path": "aiida_workgraph.engine.workgraph",
        "callable_name": "WorkGraphEngine",
        "graph_data": wgdata,
        "type": tdata["metadata"]["task_type"],
    }
    tdata["metadata"]["group_outputs"] = group_outputs
    tdata["executor"] = executor
    task = create_task(tdata)
    return task


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


def generate_tdata(
    func: Callable,
    identifier: str,
    inputs: List[Tuple[str, str]],
    outputs: List[Tuple[str, str]],
    properties: List[Tuple[str, str]],
    catalog: str,
    task_type: str,
    group_inputs: List[Tuple[str, str]] = None,
    group_outputs: List[Tuple[str, str]] = None,
    additional_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate task data for creating a task."""
    from node_graph.decorator import generate_input_sockets

    task_inputs = generate_input_sockets(
        func, inputs, properties, type_mapping=type_mapping
    )
    for input in task_inputs:
        input.setdefault("metadata", {})
        input["metadata"]["is_function_input"] = True
    task_outputs = outputs
    # add built-in sockets
    for output in builtin_outputs:
        task_outputs.append(output.copy())
    for input in builtin_inputs:
        task_inputs.append(input.copy())
    tdata = {
        "identifier": identifier,
        "metadata": {
            "task_type": task_type,
            "catalog": catalog,
            "node_class": {
                "module_path": "aiida_workgraph.task",
                "callable_name": "Task",
            },
            "group_inputs": group_inputs or [],
            "group_outputs": group_outputs or [],
        },
        "properties": properties,
        "inputs": task_inputs,
        "outputs": task_outputs,
    }
    tdata["executor"] = build_callable(func)
    if additional_data:
        tdata.update(additional_data)
    return tdata


class TaskDecoratorCollection:
    """Collection of task decorators."""

    # decorator with arguments indentifier, properties, inputs, outputs, executor
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

        if inputs:
            inputs = validate_task_inout(inputs, "inputs")

        if outputs:
            outputs = validate_task_inout(outputs, "outputs")

        def decorator(func):
            nonlocal identifier, task_type

            if identifier is None:
                identifier = func.__name__

            # Determine task_type based on AiiDA's node classes
            task_type = task_type
            if hasattr(func, "node_class"):
                task_type = task_types.get(func.node_class, task_type)
            task_outputs = outputs or [
                {"identifier": "workgraph.any", "name": "result"}
            ]
            for output in task_outputs:
                output.setdefault("metadata", {})
                output["metadata"]["is_function_output"] = True
            tdata = generate_tdata(
                func,
                identifier,
                inputs or [],
                task_outputs,
                properties or [],
                catalog,
                task_type,
            )
            task = create_task(tdata)
            task._error_handlers = error_handlers
            func.identifier = identifier
            func.task = func.node = task
            func.tdata = tdata
            return func

        return decorator

    # decorator with arguments indentifier, properties, inputs, outputs, executor
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

        outputs = outputs or []

        def decorator(func):
            nonlocal identifier, inputs, outputs

            if identifier is None:
                identifier = func.__name__
            # use cloudpickle to serialize function
            func.identifier = identifier

            task_outputs = [
                {"identifier": "workgraph.any", "name": output["name"]}
                for output in outputs
            ]
            # print(task_inputs, task_outputs)
            #
            task_type = "graph_builder"
            tdata = generate_tdata(
                func,
                identifier,
                inputs or [],
                task_outputs,
                properties or [],
                catalog,
                task_type,
                group_inputs=inputs,
                group_outputs=outputs,
            )

            task = create_task(tdata)
            func.task = func.node = task
            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def calcfunction(**kwargs: Any) -> Callable:
        def decorator(func):
            # First, apply the calcfunction decorator
            func_decorated = calcfunction(func)
            # Then, apply task decorator
            task_decorated = build_task_from_callable(
                func_decorated,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            identifier = kwargs.get("identifier", None)
            func_decorated.identifier = identifier if identifier else func.__name__
            func_decorated.task = func_decorated.node = task_decorated
            return func_decorated

        return decorator

    @staticmethod
    @nonfunctional_usage
    def workfunction(**kwargs: Any) -> Callable:
        def decorator(func):
            # First, apply the workfunction decorator
            func_decorated = workfunction(func)
            task_decorated = build_task_from_callable(
                func_decorated,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            identifier = kwargs.get("identifier", None)
            func_decorated.identifier = identifier if identifier else func.__name__
            func_decorated.task = func_decorated.node = task_decorated

            return func_decorated

        return decorator

    @staticmethod
    @nonfunctional_usage
    def pythonjob(**kwargs: Any) -> Callable:
        def decorator(func):
            # first create a task from the function
            task_decorated = build_task_from_callable(
                func,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            # then build a PythonJob task from the function task
            task_decorated, _ = build_pythonjob_task(func)
            func.identifier = "PythonJob"
            func.task = func.node = task_decorated
            task_decorated._error_handlers = kwargs.get("error_handlers", [])

            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def awaitable(**kwargs: Any) -> Callable:
        def decorator(func):
            # Then, apply task decorator
            task_decorated = build_task_from_callable(
                func,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            task_decorated.node_type = "awaitable"
            func.identifier = "awaitable"
            func.task = func.node = task_decorated
            return func

        return decorator

    @staticmethod
    @nonfunctional_usage
    def monitor(**kwargs: Any) -> Callable:
        def decorator(func):
            # Then, apply task decorator
            task_decorated = build_task_from_callable(
                func,
                inputs=kwargs.get(
                    "inputs", [{"identifier": "workgraph.any", "name": "interval"}]
                ),
                outputs=kwargs.get("outputs", []),
            )
            task_decorated.node_type = "monitor"
            func.identifier = "monitor"
            func.task = func.node = task_decorated
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
