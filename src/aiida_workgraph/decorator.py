from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida_workgraph.utils import get_executor
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida import orm
from aiida.orm.nodes.process.calculation.calcfunction import CalcFunctionNode
from aiida.orm.nodes.process.workflow.workfunction import WorkFunctionNode
from aiida.engine.processes.ports import PortNamespace
from aiida_workgraph.task import Task
from aiida_workgraph.utils import build_callable, validate_task_inout
import inspect

task_types = {
    CalcFunctionNode: "CALCFUNCTION",
    WorkFunctionNode: "WORKFUNCTION",
    CalcJob: "CALCJOB",
    WorkChain: "WORKCHAIN",
}

type_mapping = {
    "default": "workgraph.any",
    "namespace": "workgraph.namespace",
    int: "workgraph.int",
    float: "workgraph.float",
    str: "workgraph.string",
    bool: "workgraph.bool",
    orm.Int: "workgraph.aiida_int",
    orm.Float: "workgraph.aiida_float",
    orm.Str: "workgraph.aiida_string",
    orm.Bool: "workgraph.aiida_bool",
    orm.StructureData: "workgraph.aiida_structuredata",
}


def create_task(tdata):
    """Wrap create_node from node_graph to create a Task."""
    from node_graph.decorator import create_node

    tdata["type_mapping"] = type_mapping
    tdata["metadata"]["node_type"] = tdata["metadata"].pop("task_type")
    return create_node(tdata)


def add_input_recursive(
    inputs: List[List[Union[str, Dict[str, Any]]]],
    port: PortNamespace,
    args: List,
    kwargs: List,
    prefix: Optional[str] = None,
    required: bool = True,
) -> List[List[Union[str, Dict[str, Any]]]]:
    """Add input recursively."""
    if prefix is None:
        port_name = port.name
    else:
        port_name = f"{prefix}.{port.name}"
    required = port.required and required
    input_names = [input["name"] for input in inputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port_name not in input_names:
            inputs.append(
                {
                    "identifier": "workgraph.namespace",
                    "name": port_name,
                    "property": {"identifier": "workgraph.any", "default": None},
                }
            )
        if required:
            args.append(port_name)
        else:
            kwargs.append(port_name)
        for value in port.values():
            add_input_recursive(
                inputs, value, args, kwargs, prefix=port_name, required=required
            )
    else:
        if port_name not in input_names:
            # port.valid_type can be a single type or a tuple of types,
            # we only support single type for now
            if isinstance(port.valid_type, tuple) and len(port.valid_type) > 1:
                socket_type = "workgraph.any"
            if isinstance(port.valid_type, tuple) and len(port.valid_type) == 1:
                socket_type = type_mapping.get(port.valid_type[0], "workgraph.any")
            else:
                socket_type = type_mapping.get(port.valid_type, "workgraph.any")
            inputs.append({"identifier": socket_type, "name": port_name})
        if required:
            args.append(port_name)
        else:
            kwargs.append(port_name)
    return inputs


def add_output_recursive(
    outputs: List[List[Union[str, Dict[str, Any]]]],
    port: PortNamespace,
    prefix: Optional[str] = None,
    required: bool = True,
) -> List[List[Union[str, Dict[str, Any]]]]:
    """Add output recursively."""
    if prefix is None:
        port_name = port.name
    else:
        port_name = f"{prefix}.{port.name}"
    required = port.required and required
    output_names = [output["name"] for output in outputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port_name not in output_names:
            outputs.append({"identifier": "workgraph.namespace", "name": port_name})
        for value in port.values():
            add_output_recursive(outputs, value, prefix=port_name, required=required)
    else:
        if port_name not in output_names:
            outputs.append({"identifier": "workgraph.any", "name": port_name})
    return outputs


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
        executor, _ = get_executor({"module": module, "name": executor_name})
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
            tdata["executor"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        else:
            tdata["metadata"]["task_type"] = "NORMAL"
            tdata["executor"] = executor
            return build_task_from_function(executor, inputs=inputs, outputs=outputs)
    else:
        if issubclass(executor, CalcJob):
            tdata["metadata"]["task_type"] = "CALCJOB"
            tdata["executor"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        elif issubclass(executor, WorkChain):
            tdata["metadata"]["task_type"] = "WORKCHAIN"
            tdata["executor"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
    raise ValueError("The executor is not supported.")


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
) -> Task:
    """Register a task from a AiiDA component.
    For example: CalcJob, WorkChain, CalcFunction, WorkFunction."""

    tdata.setdefault("metadata", {})
    inputs = [] if inputs is None else inputs
    outputs = [] if outputs is None else outputs
    executor = tdata["executor"]
    spec = executor.spec()
    args = []
    kwargs = []
    user_defined_input_names = [input["name"] for input in inputs]
    for _key, port in spec.inputs.ports.items():
        add_input_recursive(inputs, port, args, kwargs, required=port.required)
    for _key, port in spec.outputs.ports.items():
        add_output_recursive(outputs, port, required=port.required)
    # Only check this for calcfunction and workfunction
    if inspect.isfunction(executor) and spec.inputs.dynamic:
        if hasattr(executor.process_class, "_varargs"):
            name = executor.process_class._varargs
        else:
            name = (
                executor.process_class._var_keyword
                or executor.process_class._var_positional
            )
        tdata["var_kwargs"] = name
        # if user already defined the var_args in the inputs, skip it
        if name not in [input["name"] for input in inputs]:
            inputs.append(
                {
                    "identifier": "workgraph.any",
                    "name": name,
                    "property": {"identifier": "workgraph.any", "default": {}},
                }
            )
        # When the input is dyanmic, if user defines some input names does not included in the args and kwargs,
        # which means the user define the input names manually, we must add them to the kwargs
        for key in user_defined_input_names:
            if key not in args and key not in kwargs:
                if key == name:
                    continue
                kwargs.append(key)

    # TODO In order to reload the WorkGraph from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    tdata["identifier"] = tdata.pop("identifier", tdata["executor"].__name__)
    tdata["executor"] = build_callable(executor)
    tdata["executor"]["type"] = tdata["metadata"]["task_type"]
    if tdata["metadata"]["task_type"].upper() in ["CALCFUNCTION", "WORKFUNCTION"]:
        outputs = (
            [{"identifier": "workgraph.any", "name": "result"}]
            if not outputs
            else outputs
        )
    # add built-in sockets
    outputs.append({"identifier": "workgraph.any", "name": "_outputs"})
    outputs.append({"identifier": "workgraph.any", "name": "_wait"})
    inputs.append({"identifier": "workgraph.any", "name": "_wait", "link_limit": 1e6})
    tdata["metadata"]["node_class"] = {"module": "aiida_workgraph.task", "name": "Task"}
    tdata["args"] = args
    tdata["kwargs"] = kwargs
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
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
        "executor": PythonJob,
    }
    _, tdata_py = build_task_from_AiiDA(tdata)
    tdata = deepcopy(func.tdata)
    function_inputs = [
        name for name in tdata["inputs"] if name not in ["_wait", "_outputs"]
    ]
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
    # append the kwargs of the PythonJob task to the function task
    kwargs = tdata["kwargs"]
    kwargs.extend(["computer", "command_info"])
    kwargs.extend(tdata_py["kwargs"])
    tdata["kwargs"] = kwargs
    tdata["metadata"]["task_type"] = "PYTHONJOB"
    tdata["identifier"] = "workgraph.pythonjob"
    tdata["metadata"]["node_class"] = {
        "module": "aiida_workgraph.tasks.pythonjob",
        "name": "PythonJob",
    }
    task = create_task(tdata)
    task.is_aiida_component = True
    task.function_inputs = function_inputs
    return task, tdata


def build_shelljob_task(
    nodes: dict = None, outputs: list = None, parser_outputs: list = None
) -> Task:
    """Build ShellJob with custom inputs and outputs."""
    from aiida_shell.calculations.shell import ShellJob
    from aiida_shell.parsers.shell import ShellParser
    from node_graph.socket import NodeSocket

    tdata = {
        "metadata": {"task_type": "SHELLJOB"},
        "executor": ShellJob,
    }
    _, tdata = build_task_from_AiiDA(tdata)
    # create input sockets for the nodes, if it is linked other sockets
    links = {}
    inputs = []
    nodes = {} if nodes is None else nodes
    keys = list(nodes.keys())
    for key in keys:
        inputs.append({"identifier": "workgraph.any", "name": f"nodes.{key}"})
        # input is a output of another task, we make a link
        if isinstance(nodes[key], NodeSocket):
            links[f"nodes.{key}"] = nodes[key]
            # Output socket itself is not a value, so we remove the key from the nodes
            nodes.pop(key)
    for input in inputs:
        if input["name"] not in tdata["inputs"]:
            input["list_index"] = len(tdata["inputs"]) + 1
            tdata["inputs"][input["name"]] = input
            tdata["kwargs"].append(input["name"])
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
    tdata["kwargs"].extend(["command", "resolve_command"])
    tdata["metadata"]["task_type"] = "SHELLJOB"
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata, links


def build_task_from_workgraph(wg: any) -> Task:
    """Build task from workgraph."""
    from aiida.orm.utils.serialize import serialize

    tdata = {"metadata": {"task_type": "workgraph"}}
    inputs = []
    outputs = []
    group_outputs = []
    # add all the inputs/outputs from the tasks in the workgraph
    for task in wg.tasks:
        # inputs
        inputs.append(
            {
                "identifier": "workgraph.namespace",
                "name": f"{task.name}",
            }
        )
        for socket in task.inputs:
            if socket.name == "_wait":
                continue
            inputs.append(
                {"identifier": socket.identifier, "name": f"{task.name}.{socket.name}"}
            )
        # outputs
        outputs.append(
            {
                "identifier": "workgraph.namespace",
                "name": f"{task.name}",
            }
        )
        for socket in task.outputs:
            if socket.name in ["_wait", "_outputs"]:
                continue
            outputs.append(
                {"identifier": socket.identifier, "name": f"{task.name}.{socket.name}"}
            )
            group_outputs.append(
                {
                    "name": f"{task.name}.{socket.name}",
                    "from": f"{task.name}.{socket.name}",
                }
            )
    kwargs = [input["name"] for input in inputs]
    # add built-in sockets
    outputs.append({"identifier": "workgraph.any", "name": "_outputs"})
    outputs.append({"identifier": "workgraph.any", "name": "_wait"})
    inputs.append({"identifier": "workgraph.any", "name": "_wait", "link_limit": 1e6})
    tdata["metadata"]["node_class"] = {"module": "aiida_workgraph.task", "name": "Task"}
    tdata["kwargs"] = kwargs
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    tdata["identifier"] = wg.name
    executor = {
        "module": "aiida_workgraph.engine.workgraph",
        "name": "WorkGraphEngine",
        "wgdata": serialize(wg.to_dict(store_nodes=True)),
        "type": tdata["metadata"]["task_type"],
        "is_pickle": False,
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

    args, kwargs, var_args, var_kwargs, _inputs = generate_input_sockets(
        func, inputs, properties, type_mapping=type_mapping
    )
    task_outputs = outputs
    # add built-in sockets
    _inputs.append({"identifier": "workgraph.any", "name": "_wait", "link_limit": 1e6})
    task_outputs.append({"identifier": "workgraph.any", "name": "_wait"})
    task_outputs.append({"identifier": "workgraph.any", "name": "_outputs"})
    tdata = {
        "identifier": identifier,
        "args": args,
        "kwargs": kwargs,
        "var_args": var_args,
        "var_kwargs": var_kwargs,
        "metadata": {
            "task_type": task_type,
            "catalog": catalog,
            "node_class": {"module": "aiida_workgraph.task", "name": "Task"},
            "group_inputs": group_inputs or [],
            "group_outputs": group_outputs or [],
        },
        "properties": properties,
        "inputs": _inputs,
        "outputs": task_outputs,
    }
    tdata["executor"] = build_callable(func)
    if additional_data:
        tdata.update(additional_data)
    return tdata


class TaskDecoratorCollection:
    """Collection of task decorators."""

    # decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
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
            args (list): task args
            kwargs (dict): task kwargs
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
            tdata = generate_tdata(
                func,
                identifier,
                inputs or [],
                outputs or [{"identifier": "workgraph.any", "name": "result"}],
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

    # decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
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
