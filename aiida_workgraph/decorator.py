from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from aiida_workgraph.utils import get_executor
from aiida.engine import calcfunction, workfunction, CalcJob, WorkChain
from aiida import orm
from aiida.orm.nodes.process.calculation.calcfunction import CalcFunctionNode
from aiida.orm.nodes.process.workflow.workfunction import WorkFunctionNode
from aiida.engine.processes.ports import PortNamespace
from node_graph.decorator import create_node
import cloudpickle as pickle
from aiida_workgraph.node import Node

node_types = {
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
        # so if you change the value of one port, the value of all the ports of other nodes will be changed
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
        # so if you change the value of one port, the value of all the ports of other nodes will be changed
        # consider to use None as default value
        if port_name not in output_names:
            outputs.append(["General", port_name])
        for value in port.values():
            add_output_recursive(outputs, value, prefix=port_name, required=required)
    else:
        if port_name not in output_names:
            outputs.append(["General", port_name])
    return outputs


def build_node(
    executor: Union[Callable, str],
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Node:
    """Build node from executor."""
    from aiida_workgraph.workgraph import WorkGraph

    if isinstance(executor, WorkGraph):
        return build_node_from_workgraph(executor)
    elif isinstance(executor, str):
        (
            path,
            executor_name,
        ) = executor.rsplit(".", 1)
        executor, _ = get_executor({"path": path, "name": executor_name})
    if callable(executor):
        return build_node_from_callable(executor, inputs=inputs, outputs=outputs)


def build_node_from_callable(
    executor: Callable,
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Node:
    """Build node from a callable object.
    First, check if the executor is already a node.
    If not, check if it is a function or a class.
    If it is a function, build node from function.
    If it is a class, it only supports CalcJob and WorkChain.
    """
    import inspect
    from aiida_workgraph.node import Node

    # if it is already a node, return it
    if (
        hasattr(executor, "node")
        and inspect.isclass(executor.node)
        and issubclass(executor.node, Node)
        or inspect.isclass(executor)
        and issubclass(executor, Node)
    ):
        return executor
    ndata = {}
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, "node_class", False):
            ndata["node_type"] = node_types.get(executor.node_class, "NORMAL")
            ndata["executor"] = executor
            return build_node_from_AiiDA(ndata, inputs=inputs, outputs=outputs)[0]
        else:
            ndata["node_type"] = "NORMAL"
            ndata["executor"] = executor
            return build_node_from_function(executor, inputs=inputs, outputs=outputs)
    else:
        if issubclass(executor, CalcJob):
            ndata["node_type"] = "CALCJOB"
            ndata["executor"] = executor
            return build_node_from_AiiDA(ndata, inputs=inputs, outputs=outputs)[0]
        elif issubclass(executor, WorkChain):
            ndata["node_type"] = "WORKCHAIN"
            ndata["executor"] = executor
            return build_node_from_AiiDA(ndata, inputs=inputs, outputs=outputs)[0]
    raise ValueError("The executor is not supported.")


def build_node_from_function(
    executor: Callable,
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Node:
    """Build node from function."""
    return NodeDecoratorCollection.decorator_node(inputs=inputs, outputs=outputs)(
        executor
    ).node


def build_node_from_AiiDA(
    ndata: Dict[str, Any],
    inputs: Optional[List[str]] = None,
    outputs: Optional[List[str]] = None,
) -> Node:
    """Register a node from a AiiDA component.
    For example: CalcJob, WorkChain, CalcFunction, WorkFunction."""
    from aiida_workgraph.node import Node

    # print(executor)
    inputs = [] if inputs is None else inputs
    outputs = [] if outputs is None else outputs
    executor = ndata["executor"]
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
        ndata["var_kwargs"] = name
        inputs.append(["General", name, {"property": ["General", {"default": {}}]}])
    if ndata["node_type"].upper() in ["CALCFUNCTION", "WORKFUNCTION"]:
        outputs = [["General", "result"]] if not outputs else outputs
    # print("kwargs: ", kwargs)
    # add built-in sockets
    outputs.append(["General", "_outputs"])
    outputs.append(["General", "_wait"])
    inputs.append(["General", "_wait", {"link_limit": 1e6}])
    ndata["node_class"] = Node
    ndata["args"] = args
    ndata["kwargs"] = kwargs
    ndata["inputs"] = inputs
    ndata["outputs"] = outputs
    ndata["identifier"] = ndata.pop("identifier", ndata["executor"].__name__)
    # TODO In order to reload the WorkGraph from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    executor = {
        "executor": pickle.dumps(executor),
        "type": ndata["node_type"],
        "is_pickle": True,
    }
    ndata["executor"] = executor
    node = create_node(ndata)
    return node, ndata


def build_PythonJob_node(func: Callable) -> Node:
    """Build PythonJob node from function."""
    from aiida_workgraph.calculations.python import PythonJob

    ndata = {"executor": PythonJob, "node_type": "CALCJOB"}
    _, ndata_py = build_node_from_AiiDA(ndata)
    ndata = func.ndata
    # merge the inputs and outputs from the PythonJob node to the function node
    # skip the already existed inputs and outputs
    inputs = ndata["inputs"]
    outputs = ndata["outputs"]
    for input in ndata_py["inputs"]:
        if input not in inputs:
            inputs.append(input)
    for output in ndata_py["outputs"]:
        if output not in outputs:
            outputs.append(output)
    # append the kwargs of the PythonJob node to the function node
    kwargs = ndata["kwargs"]
    kwargs.extend(ndata_py["kwargs"])
    ndata["inputs"] = inputs
    ndata["outputs"] = outputs
    ndata["kwargs"] = kwargs
    ndata["node_type"] = "PYTHONJOB"
    return create_node(ndata)


def build_node_from_workgraph(wg: any) -> Node:
    """Build node from workgraph."""
    from aiida_workgraph.node import Node

    ndata = {"node_type": "workgraph"}
    inputs = []
    outputs = []
    group_outputs = []
    # add all the inputs/outputs from the nodes in the workgraph
    for node in wg.nodes:
        # inputs
        inputs.append(["General", f"{node.name}"])
        for socket in node.inputs:
            if socket.name == "_wait":
                continue
            inputs.append(["General", f"{node.name}.{socket.name}"])
        # outputs
        outputs.append(["General", f"{node.name}"])
        for socket in node.outputs:
            if socket.name in ["_wait", "_outputs"]:
                continue
            outputs.append(["General", f"{node.name}.{socket.name}"])
            group_outputs.append(
                [f"{node.name}.{socket.name}", f"{node.name}.{socket.name}"]
            )
    kwargs = [input[1] for input in inputs]
    # add built-in sockets
    outputs.append(["General", "_outputs"])
    outputs.append(["General", "_wait"])
    inputs.append(["General", "_wait", {"link_limit": 1e6}])
    inputs.append(["General", "_code"])
    ndata["node_class"] = Node
    ndata["kwargs"] = kwargs
    ndata["inputs"] = inputs
    ndata["outputs"] = outputs
    ndata["identifier"] = wg.name
    # TODO In order to reload the WorkGraph from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    executor = {
        "executor": None,
        "wgdata": wg.to_dict(),
        "type": ndata["node_type"],
        "is_pickle": True,
    }
    ndata["executor"] = executor
    node = create_node(ndata)
    node.group_outputs = group_outputs
    return node


def serialize_function(func: Callable) -> Dict[str, Any]:
    """Serialize a function for storage or transmission."""
    import cloudpickle as pickle
    import inspect
    import textwrap

    source_code = inspect.getsource(func)
    source_code_lines = source_code.split("\n")
    function_source_code = "\n".join(source_code_lines[1:])
    function_source_code = textwrap.dedent(function_source_code)

    return {
        "executor": pickle.dumps(func),
        "type": "function",
        "is_pickle": True,
        "function_source_code": function_source_code,
    }


def generate_ndata(
    func: Callable,
    identifier: str,
    inputs: List[Tuple[str, str]],
    outputs: List[Tuple[str, str]],
    properties: List[Tuple[str, str]],
    catalog: str,
    node_type: str,
    additional_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate node data for creating a node."""
    from node_graph.decorator import generate_input_sockets
    from aiida_workgraph.node import Node

    args, kwargs, var_args, var_kwargs, _inputs = generate_input_sockets(
        func, inputs, properties
    )
    node_outputs = outputs
    # add built-in sockets
    _inputs.append(["General", "_wait", {"link_limit": 1e6}])
    _inputs.append(["General", "_code"])
    kwargs.append("_code")
    node_outputs.append(["General", "_wait"])
    node_outputs.append(["General", "_outputs"])
    ndata = {
        "node_class": Node,
        "identifier": identifier,
        "args": args,
        "kwargs": kwargs,
        "var_args": var_args,
        "var_kwargs": var_kwargs,
        "node_type": node_type,
        "properties": properties,
        "inputs": _inputs,
        "outputs": node_outputs,
        "executor": serialize_function(func),
        "catalog": catalog,
    }
    if additional_data:
        ndata.update(additional_data)
    return ndata


class NodeDecoratorCollection:
    """Collection of node decorators."""

    # decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
    @staticmethod
    def decorator_node(
        identifier: Optional[str] = None,
        node_type: str = "Normal",
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Tuple[str, str]]] = None,
        outputs: Optional[List[Tuple[str, str]]] = None,
        catalog: str = "Others",
    ) -> Callable:
        """Generate a decorator that register a function as a node.

        Attributes:
            indentifier (str): node identifier
            catalog (str): node catalog
            args (list): node args
            kwargs (dict): node kwargs
            properties (list): node properties
            inputs (list): node inputs
            outputs (list): node outputs
        """

        def decorator(func):
            nonlocal identifier, node_type

            if identifier is None:
                identifier = func.__name__

            # Determine node_type based on AiiDA's node classes
            node_type = node_type
            if hasattr(func, "node_class"):
                node_type = node_types.get(func.node_class, node_type)
            ndata = generate_ndata(
                func,
                identifier,
                inputs or [],
                outputs or [["General", "result"]],
                properties or [],
                catalog,
                node_type,
            )
            node = create_node(ndata)
            func.identifier = identifier
            func.node = node
            func.ndata = ndata
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
        """Generate a decorator that register a node group as a node.
        Attributes:
            indentifier (str): node identifier
            catalog (str): node catalog
            properties (list): node properties
            inputs (list): node inputs
            outputs (list): node outputs
        """

        outputs = outputs or []

        def decorator(func):
            nonlocal identifier, inputs, outputs

            if identifier is None:
                identifier = func.__name__
            # use cloudpickle to serialize function
            func.identifier = identifier

            node_outputs = [["General", output[1]] for output in outputs]
            # print(node_inputs, node_outputs)
            #
            node_type = "graph_builder"
            ndata = generate_ndata(
                func,
                identifier,
                inputs or [],
                node_outputs,
                properties or [],
                catalog,
                node_type,
            )
            node = create_node(ndata)
            node.group_inputs = inputs
            node.group_outputs = outputs
            func.node = node
            return func

        return decorator

    @staticmethod
    def calcfunction(**kwargs: Any) -> Callable:
        def decorator(func):
            # First, apply the calcfunction decorator
            func_decorated = calcfunction(func)
            # Then, apply node decorator
            node_decorated = build_node_from_callable(
                func_decorated,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            identifier = kwargs.get("identifier", None)
            func_decorated.identifier = identifier if identifier else func.__name__
            func_decorated.node = node_decorated
            return func_decorated

        return decorator

    @staticmethod
    def workfunction(**kwargs: Any) -> Callable:
        def decorator(func):
            # First, apply the workfunction decorator
            func_decorated = workfunction(func)
            node_decorated = build_node_from_callable(
                func_decorated,
                inputs=kwargs.get("inputs", []),
                outputs=kwargs.get("outputs", []),
            )
            identifier = kwargs.get("identifier", None)
            func_decorated.identifier = identifier if identifier else func.__name__
            func_decorated.node = node_decorated

            return func_decorated

        return decorator

    # Making decorator_node accessible as 'node'
    node = decorator_node

    # Making decorator_graph_builder accessible as 'graph_builder'
    graph_builder = decorator_graph_builder

    def __call__(self, *args, **kwargs):
        # This allows using '@node' to directly apply the decorator_node functionality
        return self.decorator_node(*args, **kwargs)


node = NodeDecoratorCollection()
