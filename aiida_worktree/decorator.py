from typing import Any
from aiida_worktree.utils import get_executor
from aiida.engine.processes.functions import calcfunction, workfunction
from aiida.engine.processes.calcjobs import CalcJob
from aiida.engine.processes.workchains import WorkChain
from aiida.orm.nodes.process.calculation.calcfunction import CalcFunctionNode
from aiida.orm.nodes.process.workflow.workfunction import WorkFunctionNode
from aiida.engine.processes.ports import PortNamespace


def add_input_recursive(inputs, port, prefix=None):
    """Add input recursively."""
    if prefix is None:
        port_name = port.name
    else:
        port_name = f"{prefix}.{port.name}"
    if isinstance(port, PortNamespace):
        # TODO the default value is {} could cause problem, because the address of the dict is the same,
        # so if you change the value of one port, the value of all the ports of other nodes will be changed
        # consider to use None as default value
        inputs.append(
            ["General", port_name, {"property": ["General", {"default": {}}]}]
        )
        for value in port.values():
            add_input_recursive(inputs, value, prefix=port_name)
    else:
        inputs.append(["General", port_name])
    return inputs


def build_node(ndata):
    """Register a node from a AiiDA component."""
    from node_graph.decorator import create_node
    from aiida_worktree.node import Node
    import cloudpickle as pickle
    import inspect

    path, executor_name, = ndata.pop(
        "path"
    ).rsplit(".", 1)
    ndata["executor"] = {"path": path, "name": executor_name}
    executor, type = get_executor(ndata["executor"])
    # print(executor)
    if inspect.isfunction(executor):
        # calcfunction and workfunction
        if getattr(executor, "node_class", False):
            if executor.node_class is CalcFunctionNode:
                ndata["node_type"] = "calcfunction"
            elif executor.node_class is WorkFunctionNode:
                ndata["node_type"] = "workfunction"
        else:
            ndata["node_type"] = "normal"
    else:
        if issubclass(executor, CalcJob):
            ndata["node_type"] = "calcjob"
        elif issubclass(executor, WorkChain):
            ndata["node_type"] = "workchain"
        else:
            ndata["node_type"] = "normal"
    inputs = []
    outputs = []
    spec = executor.spec()
    for _key, port in spec.inputs.ports.items():
        add_input_recursive(inputs, port)
    kwargs = [input[1] for input in inputs]
    for _key, port in spec.outputs.ports.items():
        outputs.append(["General", port.name])
    if spec.inputs.dynamic:
        ndata["var_kwargs"] = spec.inputs.dynamic
        inputs.append(
            ["General", spec.inputs.dynamic, {"property": ["General", {"default": {}}]}]
        )
    # print("kwargs: ", kwargs)
    ndata["node_class"] = Node
    ndata["kwargs"] = kwargs
    ndata["inputs"] = inputs
    ndata["outputs"] = outputs
    ndata["identifier"] = ndata.pop("identifier", ndata["executor"]["name"])
    # TODO In order to reload the WorkTree from process, "is_pickle" should be True
    # so I pickled the function here, but this is not necessary
    # we need to update the node_graph to support the path and name of the function
    executor = {
        "executor": pickle.dumps(executor),
        "type": "function",
        "is_pickle": True,
    }
    ndata["executor"] = executor
    node = create_node(ndata)
    return node


# decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
def decorator_node(
    identifier=None,
    node_type="Normal",
    properties=None,
    inputs=None,
    outputs=None,
    catalog="Others",
    executor_type="function",
):
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
    from aiida_worktree.node import Node

    properties = properties or []
    inputs = inputs or []
    outputs = outputs or [["General", "result"]]

    def decorator(func):
        import cloudpickle as pickle
        from node_graph.decorator import generate_input_sockets, create_node

        nonlocal identifier

        if identifier is None:
            identifier = func.__name__
        # use cloudpickle to serialize function
        executor = {
            "executor": pickle.dumps(func),
            "type": "function",
            "is_pickle": True,
        }
        # Get the args and kwargs of the function
        args, kwargs, var_args, var_kwargs, _inputs = generate_input_sockets(
            func, inputs, properties
        )
        # I don't know why isinstance(func.node_class, CalcFunctionNode) is False
        # print("func: ", func)
        # print("node_class: ", func.node_class)
        if hasattr(func, "node_class"):
            if func.node_class is CalcFunctionNode:
                node_type = "calcfunction"
            elif func.node_class is WorkFunctionNode:
                node_type = "workfunction"
        else:
            node_type = "Normal"
        ndata = {
            "node_class": Node,
            "identifier": identifier,
            "node_type": node_type,
            "args": args,
            "kwargs": kwargs,
            "var_args": var_args,
            "var_kwargs": var_kwargs,
            "properties": properties,
            "inputs": _inputs,
            "outputs": outputs,
            "executor": executor,
            "catalog": catalog,
        }
        node = create_node(ndata)
        func.identifier = identifier
        func.node = node
        return func

    return decorator


# decorator with arguments indentifier, args, kwargs, properties, inputs, outputs, executor
def decorator_node_group(
    identifier=None,
    properties=None,
    inputs=None,
    outputs=None,
    catalog="Others",
    executor_type="function",
):
    """Generate a decorator that register a node group as a node.

    Attributes:
        indentifier (str): node identifier
        catalog (str): node catalog
        properties (list): node properties
        inputs (list): node inputs
        outputs (list): node outputs
    """
    from aiida_worktree.node import Node

    properties = properties or []
    inputs = inputs or []
    outputs = outputs or []

    def decorator(func):
        import cloudpickle as pickle
        from node_graph.decorator import generate_input_sockets, create_node

        nonlocal identifier, inputs, outputs

        if identifier is None:
            identifier = func.__name__
        # use cloudpickle to serialize function
        func.identifier = identifier
        func.group_outputs = outputs
        executor = {
            "executor": pickle.dumps(func),
            "type": executor_type,
            "is_pickle": True,
        }
        # Get the args and kwargs of the function
        args, kwargs, var_args, var_kwargs, _inputs = generate_input_sockets(
            func, inputs, properties
        )
        # nt = func()
        # inputs = [[nt.nodes[input[0]].inputs[input[1]].identifier, input[2]] for input in group_inputs]
        # outputs = [[nt.nodes[output[0]].outputs[output[1]].identifier, output[2]] for output in group_outputs]
        # node_inputs = [["General", input[2]] for input in inputs]
        node_outputs = [["General", output[1]] for output in outputs]
        # print(node_inputs, node_outputs)
        #
        node_type = "worktree"
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
            "executor": executor,
            "catalog": catalog,
        }
        node = create_node(ndata)
        node.group_inputs = inputs
        node.group_outputs = outputs
        func.node = node
        return func

    return decorator


class NodeDecoratorCollection:
    """Collection of node decorators."""

    node = staticmethod(decorator_node)
    group = staticmethod(decorator_node_group)

    __call__: Any = node  # Alias '@node' to '@node.node'.

    @staticmethod
    def calcfunction(**kwargs):
        def decorator(func):
            # First, apply the calcfunction decorator
            calcfunc_decorated = calcfunction(func)
            # Then, apply node decorator
            node_decorated = node(**kwargs)(calcfunc_decorated)

            return node_decorated

        return decorator

    @staticmethod
    def workfunction(**kwargs):
        def decorator(func):
            # First, apply the workfunction decorator
            calcfunc_decorated = workfunction(func)
            node_decorated = node(**kwargs)(calcfunc_decorated)

            return node_decorated

        return decorator


node = NodeDecoratorCollection()
