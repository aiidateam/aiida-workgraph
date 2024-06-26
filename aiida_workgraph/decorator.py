from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida_workgraph.utils import get_executor
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida import orm
from aiida.orm.nodes.process.calculation.calcfunction import CalcFunctionNode
from aiida.orm.nodes.process.workflow.workfunction import WorkFunctionNode
from aiida.engine.processes.ports import PortNamespace
import cloudpickle as pickle
from aiida_workgraph.task import Task

task_types = {
    CalcFunctionNode: "CALCFUNCTION",
    WorkFunctionNode: "WORKFUNCTION",
    CalcJob: "CALCJOB",
    WorkChain: "WORKCHAIN",
}

aiida_socket_maping = {
    orm.Int: "AiiDAInt",
    orm.Float: "AiiDAFloat",
    orm.Str: "AiiDAString",
    orm.Bool: "AiiDABool",
}


def create_task(tdata):
    """Wrap create_node from node_graph to create a Task."""
    from node_graph.decorator import create_node

    tdata["node_type"] = tdata.pop("task_type")
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
    input_names = [input[1] for input in inputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port_name not in input_names:
            inputs.append(
                ["General", port_name, {"property": ["General", {"default": {}}]}]
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
                socket_type = "General"
            if isinstance(port.valid_type, tuple) and len(port.valid_type) == 1:
                socket_type = aiida_socket_maping.get(port.valid_type[0], "General")
            else:
                socket_type = aiida_socket_maping.get(port.valid_type, "General")
            inputs.append([socket_type, port_name])
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
    output_names = [output[1] for output in outputs]
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other tasks will be changed
        # consider to use None as default value
        if port_name not in output_names:
            outputs.append(["General", port_name])
        for value in port.values():
            add_output_recursive(outputs, value, prefix=port_name, required=required)
    else:
        if port_name not in output_names:
            outputs.append(["General", port_name])
    return outputs


def build_task(
    executor: Union[Callable, str],
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Task:
    """Build task from executor."""
    from aiida_workgraph.workgraph import WorkGraph

    if isinstance(executor, WorkGraph):
        return build_task_from_workgraph(executor)
    elif isinstance(executor, str):
        (
            path,
            executor_name,
        ) = executor.rsplit(".", 1)
        executor, _ = get_executor({"path": path, "name": executor_name})
    if callable(executor):
        return build_task_from_callable(executor, inputs=inputs, outputs=outputs)


def build_task_from_callable(
    executor: Callable,
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Task:
    """Build task from a callable object.
    First, check if the executor is already a task.
    If not, check if it is a function or a class.
    If it is a function, build task from function.
    If it is a class, it only supports CalcJob and WorkChain.
    """
    import inspect
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
    tdata = {}
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, "node_class", False):
            tdata["task_type"] = task_types.get(executor.node_class, "NORMAL")
            tdata["executor"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        else:
            tdata["task_type"] = "NORMAL"
            tdata["executor"] = executor
            return build_task_from_function(executor, inputs=inputs, outputs=outputs)
    else:
        if issubclass(executor, CalcJob):
            tdata["task_type"] = "CALCJOB"
            tdata["executor"] = executor
            return build_task_from_AiiDA(tdata, inputs=inputs, outputs=outputs)[0]
        elif issubclass(executor, WorkChain):
            tdata["task_type"] = "WORKCHAIN"
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
    from aiida_workgraph.task import Task

    # print(executor)
    inputs = [] if inputs is None else inputs
    outputs = [] if outputs is None else outputs
    executor = tdata["executor"]
    spec = executor.spec()
    args = []
    kwargs = []
    for _key, port in spec.inputs.ports.items():
        add_input_recursive(inputs, port, args, kwargs, required=port.required)
    for _key, port in spec.outputs.ports.items():
        add_output_recursive(outputs, port, required=port.required)
    if spec.inputs.dynamic:
        if hasattr(executor.process_class, "_varargs"):
            name = executor.process_class._varargs
        else:
            name = (
                executor.process_class._var_keyword
                or executor.process_class._var_positional
            )
        tdata["var_kwargs"] = name
        inputs.append(["General", name, {"property": ["General", {"default": {}}]}])
    if tdata["task_type"].upper() in ["CALCFUNCTION", "WORKFUNCTION"]:
        outputs = [["General", "result"]] if not outputs else outputs
    # print("kwargs: ", kwargs)
    # add built-in sockets
    outputs.append(["General", "_outputs"])
    outputs.append(["General", "_wait"])
    inputs.append(["General", "_wait", {"link_limit": 1e6}])
    tdata["node_class"] = Task
    tdata["args"] = args
    tdata["kwargs"] = kwargs
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    tdata["identifier"] = tdata.pop("identifier", tdata["executor"].__name__)
    # TODO In order to reload the WorkGraph from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    executor = {
        "executor": pickle.dumps(executor),
        "type": tdata["task_type"],
        "is_pickle": True,
    }
    tdata["executor"] = executor
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata


def build_python_task(func: Callable) -> Task:
    """Build PythonJob task from function."""
    from aiida_workgraph.calculations.python import PythonJob
    from copy import deepcopy

    tdata = {"executor": PythonJob, "task_type": "CALCJOB"}
    _, tdata_py = build_task_from_AiiDA(tdata)
    tdata = deepcopy(func.tdata)
    # merge the inputs and outputs from the PythonJob task to the function task
    # skip the already existed inputs and outputs
    inputs = tdata["inputs"]
    inputs.extend(
        [
            ["String", "computer"],
            ["String", "code_label"],
            ["String", "code_path"],
            ["String", "prepend_text"],
        ]
    )
    outputs = tdata["outputs"]
    for input in tdata_py["inputs"]:
        if input not in inputs:
            inputs.append(input)
    for output in tdata_py["outputs"]:
        if output not in outputs:
            outputs.append(output)
    # change "copy_files" link_limit to 1e6
    for input in inputs:
        if input[1] == "copy_files":
            input[2].update({"link_limit": 1e6})
    # append the kwargs of the PythonJob task to the function task
    kwargs = tdata["kwargs"]
    kwargs.extend(["computer", "code_label", "code_path", "prepend_text"])
    kwargs.extend(tdata_py["kwargs"])
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    tdata["kwargs"] = kwargs
    tdata["task_type"] = "PYTHONJOB"
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata


def build_shell_task(
    nodes: dict = None, outputs: list = None, parser_outputs: list = None
) -> Task:
    """Build ShellJob with custom inputs and outputs."""
    from aiida_shell.calculations.shell import ShellJob
    from aiida_shell.parsers.shell import ShellParser
    from node_graph.socket import NodeSocket

    tdata = {"executor": ShellJob, "task_type": "SHELLJOB"}
    _, tdata = build_task_from_AiiDA(tdata)
    # create input sockets for the nodes, if it is linked other sockets
    links = {}
    inputs = []
    nodes = {} if nodes is None else nodes
    keys = list(nodes.keys())
    for key in keys:
        inputs.append(["General", f"nodes.{key}"])
        # input is a output of another task, we make a link
        if isinstance(nodes[key], NodeSocket):
            links[f"nodes.{key}"] = nodes[key]
            # Output socket itself is not a value, so we remove the key from the nodes
            nodes.pop(key)
    for input in inputs:
        if input not in tdata["inputs"]:
            tdata["inputs"].append(input)
            tdata["kwargs"].append(input[1])
    # Extend the outputs
    tdata["outputs"].extend([["General", "stdout"], ["General", "stderr"]])
    outputs = [] if outputs is None else outputs
    parser_outputs = [] if parser_outputs is None else parser_outputs
    outputs = [["General", ShellParser.format_link_label(output)] for output in outputs]
    outputs.extend(parser_outputs)
    # add user defined outputs
    for output in outputs:
        if output not in tdata["outputs"]:
            tdata["outputs"].append(output)
    #
    tdata["identifier"] = "ShellJob"
    tdata["inputs"].extend(
        [
            ["General", "command"],
            ["General", "resolve_command"],
        ]
    )
    tdata["kwargs"].extend(["command", "resolve_command"])
    tdata["task_type"] = "SHELLJOB"
    task = create_task(tdata)
    task.is_aiida_component = True
    return task, tdata, links


def build_task_from_workgraph(wg: any) -> Task:
    """Build task from workgraph."""
    from aiida_workgraph.task import Task

    tdata = {"task_type": "workgraph"}
    inputs = []
    outputs = []
    group_outputs = []
    # add all the inputs/outputs from the tasks in the workgraph
    for task in wg.tasks:
        # inputs
        inputs.append(["General", f"{task.name}"])
        for socket in task.inputs:
            if socket.name == "_wait":
                continue
            inputs.append(["General", f"{task.name}.{socket.name}"])
        # outputs
        outputs.append(["General", f"{task.name}"])
        for socket in task.outputs:
            if socket.name in ["_wait", "_outputs"]:
                continue
            outputs.append(["General", f"{task.name}.{socket.name}"])
            group_outputs.append(
                [f"{task.name}.{socket.name}", f"{task.name}.{socket.name}"]
            )
    kwargs = [input[1] for input in inputs]
    # add built-in sockets
    outputs.append(["General", "_outputs"])
    outputs.append(["General", "_wait"])
    inputs.append(["General", "_wait", {"link_limit": 1e6}])
    tdata["node_class"] = Task
    tdata["kwargs"] = kwargs
    tdata["inputs"] = inputs
    tdata["outputs"] = outputs
    tdata["identifier"] = wg.name
    # TODO In order to reload the WorkGraph from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    executor = {
        "executor": None,
        "wgdata": wg.to_dict(),
        "type": tdata["task_type"],
        "is_pickle": True,
    }
    tdata["executor"] = executor
    task = create_task(tdata)
    task.group_outputs = group_outputs
    return task


def get_required_imports(func):
    """Retrieve type hints and the corresponding module"""
    from typing import get_type_hints, _SpecialForm

    type_hints = get_type_hints(func)
    imports = {}

    def add_imports(type_hint):
        if isinstance(
            type_hint, _SpecialForm
        ):  # Handle special forms like Any, Union, Optional
            module_name = "typing"
            type_name = type_hint._name or str(type_hint)
        elif hasattr(
            type_hint, "__origin__"
        ):  # This checks for higher-order types like List, Dict
            module_name = type_hint.__module__
            type_name = type_hint._name
            for arg in type_hint.__args__:
                if arg is type(None):  # noqa: E721
                    continue
                add_imports(arg)  # Recursively add imports for each argument
        elif hasattr(type_hint, "__module__"):
            module_name = type_hint.__module__
            type_name = type_hint.__name__
        else:
            return  # If no module or origin, we can't import it, e.g., for literals

        if module_name not in imports:
            imports[module_name] = set()
        imports[module_name].add(type_name)

    for _, type_hint in type_hints.items():
        add_imports(type_hint)

    return imports


def serialize_function(func: Callable) -> Dict[str, Any]:
    """Serialize a function for storage or transmission."""
    import cloudpickle as pickle
    import inspect
    import textwrap

    source_code = inspect.getsource(func)
    source_code_lines = source_code.split("\n")
    # we need save the source code explicitly, because in the case of jupyter notebook,
    # the source code is not saved in the pickle file
    function_source_code = "\n".join(source_code_lines[1:])
    function_source_code = textwrap.dedent(function_source_code)
    # we also need to include the necessary imports for the types used in the type hints.
    try:
        required_imports = get_required_imports(func)
    except Exception as e:
        required_imports = {}
        print(f"Failed to get required imports for function {func.__name__}: {e}")
    # Generate import statements
    import_statements = "\n".join(
        f"from {module} import {', '.join(types)}"
        for module, types in required_imports.items()
    )
    return {
        "executor": pickle.dumps(func),
        "type": "function",
        "is_pickle": True,
        "function_name": func.__name__,
        "function_source_code": function_source_code,
        "import_statements": import_statements,
    }


def generate_tdata(
    func: Callable,
    identifier: str,
    inputs: List[Tuple[str, str]],
    outputs: List[Tuple[str, str]],
    properties: List[Tuple[str, str]],
    catalog: str,
    task_type: str,
    additional_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate task data for creating a task."""
    from node_graph.decorator import generate_input_sockets
    from aiida_workgraph.task import Task

    args, kwargs, var_args, var_kwargs, _inputs = generate_input_sockets(
        func, inputs, properties
    )
    task_outputs = outputs
    # add built-in sockets
    _inputs.append(["General", "_wait", {"link_limit": 1e6}])
    task_outputs.append(["General", "_wait"])
    task_outputs.append(["General", "_outputs"])
    tdata = {
        "node_class": Task,
        "identifier": identifier,
        "args": args,
        "kwargs": kwargs,
        "var_args": var_args,
        "var_kwargs": var_kwargs,
        "task_type": task_type,
        "properties": properties,
        "inputs": _inputs,
        "outputs": task_outputs,
        "executor": serialize_function(func),
        "catalog": catalog,
    }
    if additional_data:
        tdata.update(additional_data)
    return tdata


class TaskDecoratorCollection:
    """Collection of task decorators."""

    # decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
    @staticmethod
    def decorator_task(
        identifier: Optional[str] = None,
        task_type: str = "Normal",
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Tuple[str, str]]] = None,
        outputs: Optional[List[Tuple[str, str]]] = None,
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
                outputs or [["General", "result"]],
                properties or [],
                catalog,
                task_type,
            )
            task = create_task(tdata)
            func.identifier = identifier
            func.task = func.node = task
            func.tdata = tdata
            return func

        return decorator

    # decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
    @staticmethod
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

            task_outputs = [["General", output[1]] for output in outputs]
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
            )
            task = create_task(tdata)
            task.group_inputs = inputs
            task.group_outputs = outputs
            func.task = func.node = task
            return func

        return decorator

    @staticmethod
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

    # Making decorator_task accessible as 'task'
    task = decorator_task

    # Making decorator_graph_builder accessible as 'graph_builder'
    graph_builder = decorator_graph_builder

    def __call__(self, *args, **kwargs):
        # This allows using '@task' to directly apply the decorator_task functionality
        return self.decorator_task(*args, **kwargs)


task = TaskDecoratorCollection()
