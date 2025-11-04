from __future__ import annotations

from typing import Any, Dict, Optional, Union, Callable, List
from aiida.engine.processes import Process
from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.engine.runners import Runner
from aiida_workgraph.config import task_types
from aiida.engine import CalcJob, WorkChain
from aiida_pythonjob import PythonJob
from aiida_pythonjob.calculations.pyfunction import PyFunction
from aiida_shell.calculations.shell import ShellJob
import inspect
import yaml
from node_graph.socket import TaggedValue
from node_graph.socket_spec import SocketSpec
from aiida.orm.utils.serialize import serialize
from aiida_workgraph.orm.utils import deserialize_safe
from copy import deepcopy


def inspect_aiida_component_type(executor: Callable) -> str:
    task_type = None
    if isinstance(executor, type):
        if executor == PythonJob:
            task_type = 'PYTHONJOB'
        elif executor == PyFunction:
            task_type = 'PYFUNCTION'
        elif executor == ShellJob:
            task_type = 'SHELLJOB'
        elif issubclass(executor, CalcJob):
            task_type = task_types[CalcJob]
        elif issubclass(executor, WorkChain):
            task_type = task_types[WorkChain]
    elif inspect.isfunction(executor):
        if getattr(executor, 'node_class', False):
            task_type = task_types[executor.node_class]
    return task_type


def get_nested_dict(d: Dict, name: str, **kwargs) -> Any:
    """Get the value from a nested dictionary.
    If default is provided, return the default value if the key is not found.
    Otherwise, raise ValueError.
    For example:
    d = {"base": {"pw": {"parameters": 2}}}
    name = "base.pw.parameters"
    """
    from aiida.orm.utils.managers import NodeLinksManager

    keys = name.split('.')
    current = d
    for key in keys:
        if key not in current:
            if 'default' in kwargs:
                return kwargs.get('default')
            else:
                if isinstance(current, dict):
                    avaiable_keys = current.keys()
                elif isinstance(current, NodeLinksManager):
                    avaiable_keys = list(current._get_keys())
                else:
                    avaiable_keys = []
                raise ValueError(f'{name} not exist. Available keys: {avaiable_keys}')
        current = current[key]
    return current


def merge_dicts(dict1: Any, dict2: Any) -> Any:
    """Recursively merges two dictionaries."""
    for key, value in dict2.items():
        if key in dict1 and isinstance(dict1[key], dict) and isinstance(value, dict):
            # Recursively merge dictionaries
            dict1[key] = merge_dicts(dict1[key], value)
        else:
            # Overwrite or add the key
            dict1[key] = value
    return dict1


def update_nested_dict(base: Optional[Dict[str, Any]], key_path: str, value: Any) -> None:
    """
    Update or create a nested dictionary structure based on a dotted key path.

    This function allows updating a nested dictionary or creating one if `d` is `None`.
    Given a dictionary and a key path (e.g., "base.pw.parameters"), it will traverse
    or create the necessary nested structure to set the provided value at the specified
    key location. If intermediate dictionaries do not exist, they will be created.
    If the resulting dictionary is empty, it is set to `None`.

    Args:
        base (Dict[str, Any] | None): The dictionary to update, which can be `None`.
                                   If `None`, an empty dictionary will be created.
        key (str): A dotted key path string representing the nested structure.
        value (Any): The value to set at the specified key.

    Example:
        base = None
        key = "scf.pw.parameters"
        value = 2
        After running:
            update_nested_dict(d, key, value)
        The result will be:
            base = {"scf": {"pw": {"parameters": 2}}}

    Edge Case:
        If the resulting dictionary is empty after the update, it will be set to `None`.

    """

    if base is None:
        base = {}
    keys = key_path.split('.')
    current_key = keys[0]
    if len(keys) == 1:
        # Base case: Merge dictionaries or set the value directly.
        if isinstance(base.get(current_key), dict) and isinstance(value, dict):
            base[current_key] = merge_dicts(base[current_key], value)
        else:
            base[current_key] = value
    else:
        # Recursive case: Ensure the key exists and is a dictionary, then recurse.
        if current_key not in base or not isinstance(base[current_key], dict):
            base[current_key] = {}
        base[current_key] = update_nested_dict(base[current_key], '.'.join(keys[1:]), value)

    return base


def update_nested_dict_with_special_keys(data: Dict[str, Any]) -> Dict[str, Any]:
    """Update the nested dictionary with special keys like "base.pw.parameters"."""
    # Remove None

    data = {k: v for k, v in data.items() if v is not None}
    #
    special_keys = [k for k in data.keys() if '.' in k]
    for key in special_keys:
        value = data.pop(key)
        update_nested_dict(data, key, value)
    return data


def generate_provenance_graph(pk: int, output: str = None, width: str = '100%', height: str = '600px') -> Any:
    """Generate the node graph for the given node pk.
    If in Jupyter, return the graphviz object.
    Otherwise, save the graph to an HTML file.
    """

    from aiida.tools.visualization import Graph
    from aiida import orm
    from IPython.display import IFrame
    import pathlib
    from .svg_to_html import svg_to_html

    in_jupyter = False
    try:
        from IPython import get_ipython

        if get_ipython() is not None:
            in_jupyter = True
    except NameError:
        pass

    graph = Graph()
    calc_node = orm.load_node(pk)
    graph.recurse_ancestors(calc_node, annotate_links='both')
    graph.recurse_descendants(calc_node, annotate_links='both')
    g = graph.graphviz
    if not in_jupyter:
        html_content = svg_to_html(g._repr_image_svg_xml(), width, height)
        if output is None:
            pathlib.Path('html').mkdir(exist_ok=True)
            output = f'html/node_graph_{pk}.html'
        with open(output, 'w') as f:
            f.write(html_content)
        return IFrame(output, width=width, height=height)
    return g


def get_dict_from_builder(builder: Any) -> Dict:
    """Transform builder to pure dict."""
    from aiida.engine.processes.builder import ProcessBuilderNamespace

    if isinstance(builder, ProcessBuilderNamespace):
        return {k: get_dict_from_builder(v) for k, v in builder.items()}
    else:
        return builder


def clean_pickled_task_executor(tdata: Dict[str, Any]) -> None:
    """Clean the pickled executor in the task data."""
    from node_graph.executor import RuntimeExecutor
    from aiida_workgraph.executors.builtins import UnavailableExecutor

    # spec
    if 'spec' in tdata:
        executor = tdata['spec'].get('executor', {})
        if executor.get('mode', '') == 'pickled_callable':
            tdata['spec']['executor'] = RuntimeExecutor.from_callable(UnavailableExecutor).to_dict()
        if executor.get('mode', '') == 'graph':
            wgdata = executor['graph_data']
            for task in wgdata['tasks'].values():
                clean_pickled_task_executor(task)
    # error handler
    for name, handler in tdata.get('error_handlers', {}).items():
        if handler.get('mode', '') == 'pickled_callable':
            tdata['error_handlers'][name] = RuntimeExecutor.from_callable(UnavailableExecutor).to_dict()


def save_workgraph_data(node: Union[int, orm.Node], inputs: Dict[str, Any]) -> None:
    from aiida_workgraph.engine.workgraph import WorkGraphSpec

    inputs = shallow_copy_nested_dict(inputs)
    wgdata = inputs.pop(WorkGraphSpec.WORKGRAPH_DATA_KEY, {})
    task_states = {}
    task_processes = {}
    task_actions = {}
    short_wgdata = workgraph_to_short_json(wgdata)
    for name, task in wgdata['tasks'].items():
        task_states[name] = task['state']
        task_processes[name] = task['process']
        task_actions[name] = task['action']
        # clean pickled executor before save to database
        clean_pickled_task_executor(task)
    node.task_states = task_states
    node.task_processes = task_processes
    node.task_actions = task_actions
    node.workgraph_data = wgdata
    node.workgraph_data_short = short_wgdata
    node.workgraph_error_handlers = wgdata.pop('error_handlers', {})
    graph_inputs = dict(inputs.pop('graph_inputs', {}))
    tasks = dict(inputs.pop('tasks', {}))
    tasks['graph_inputs'] = graph_inputs
    node.task_inputs = serialize(tasks)


def restore_workgraph_data_from_raw_inputs(raw_inputs: Dict[str, Any]) -> Dict[str, Any]:
    """Restore the workgraph data from the raw inputs."""
    from aiida_workgraph.engine.workgraph import WorkGraphSpec

    raw_inputs = dict(raw_inputs)
    wgdata = dict(raw_inputs.pop(WorkGraphSpec.WORKGRAPH_DATA_KEY, {}))
    task_inputs = dict(raw_inputs.pop('tasks', {}))
    graph_inputs = dict(raw_inputs.pop('graph_inputs', {}))
    task_inputs['graph_inputs'] = graph_inputs
    for name, data in task_inputs.items():
        wgdata['tasks'][name]['inputs'] = data
    return wgdata


def load_workgraph_data(node: Union[int, orm.Node]) -> Optional[Dict[str, Any]]:
    """
    Get the workgraph data from the given process node.
    """
    from aiida.orm import load_node
    from aiida_workgraph.engine.workgraph import WorkGraphSpec

    if isinstance(node, int):
        node = load_node(node)
    wgdata = node.base.attributes.get(WorkGraphSpec.WORKGRAPH_DATA_KEY)
    try:
        task_inputs = deserialize_safe(node.task_inputs or '')
    except (yaml.constructor.ConstructorError, yaml.YAMLError):
        print(
            'Info: could not deserialize inputs.The workgraph is still loaded and you can inspect tasks and outputs. '
        )
        task_inputs = {}

    for name, data in task_inputs.items():
        wgdata['tasks'][name]['inputs'] = data
    wgdata['error_handlers'] = node.workgraph_error_handlers
    return wgdata


def get_parent_workgraphs(pk: int) -> List[List[str, int]]:
    """Get the list of parent workgraphs.
    Use aiida incoming links to find the parent workgraphs.
    the parent workgraph is the workgraph that has a link (type CALL_WORK) to the current workgraph.
    """
    from aiida import orm
    from aiida.common.links import LinkType

    node = orm.load_node(pk)
    parent_workgraphs = [[node.process_label, node.pk]]
    links = node.base.links.get_incoming(link_type=LinkType.CALL_WORK).all()
    if len(links) > 0:
        parent_workgraphs.extend(get_parent_workgraphs(links[0].node.pk))
    return parent_workgraphs


def get_processes_latest(
    pk: int, task_name: str = None, item_type: str = 'task'
) -> Dict[str, Dict[str, Union[int, str]]]:
    """Get the latest info of all tasks from the process."""
    import aiida
    from aiida_workgraph.orm.workgraph import WorkGraphNode

    tasks = {}
    if pk is None:
        return tasks
    node = aiida.orm.load_node(pk)
    if item_type == 'called_process':
        # fetch the process that called by the workgraph
        for link in node.base.links.get_outgoing().all():
            if isinstance(link.node, aiida.orm.ProcessNode):
                tasks[f'{link.link_label}-{link.node.pk}'] = {
                    'pk': link.node.pk,
                    'process_type': link.node.process_type,
                    'state': link.node.process_state.value,
                    'ctime': link.node.ctime,
                    'mtime': link.node.mtime,
                }
    elif item_type == 'task':
        if not isinstance(node, WorkGraphNode):
            return tasks
        task_states = node.task_states
        task_processes = node.task_processes
        task_names = [task_name] if task_name else task_states.keys()
        for name in task_names:
            state = task_states[name]
            task_process = deserialize_safe(task_processes.get(name, ''))
            tasks[name] = {
                'pk': task_process.pk if task_process else None,
                'process_type': task_process.process_type if task_process else '',
                'state': state,
                'ctime': task_process.ctime if task_process else None,
                'mtime': task_process.mtime if task_process else None,
            }

    return tasks


def get_or_create_code(
    computer: str = 'localhost',
    code_label: str = 'python3',
    code_path: str = None,
    prepend_text: str = '',
):
    """Try to load code, create if not exit."""
    from aiida.orm.nodes.data.code.installed import InstalledCode

    try:
        return orm.load_code(f'{code_label}@{computer}')
    except NotExistent:
        description = f'Code on computer: {computer}'
        computer = orm.load_computer(computer)
        code_path = code_path or code_label
        code = InstalledCode(
            computer=computer,
            label=code_label,
            description=description,
            filepath_executable=code_path,
            default_calc_job_plugin='workgraph.python',
            prepend_text=prepend_text,
        )

        code.store()
        return code


def create_and_pause_process(
    runner: Runner = None,
    process_class: Callable = None,
    inputs: dict = None,
    state_msg: str = '',
) -> Process:
    from aiida.engine.utils import instantiate_process

    process_inited = instantiate_process(runner, process_class, **inputs)
    process_inited.pause(state_msg)
    process_inited.runner.persister.save_checkpoint(process_inited)
    process_inited.close()
    runner.controller.continue_process(process_inited.pid, nowait=True, no_reply=True)
    return process_inited


def get_raw_value(identifier, value: Any) -> Any:
    """Get the raw value from a Data node."""
    if identifier in [
        'workgraph.int',
        'workgraph.float',
        'workgraph.string',
        'workgraph.bool',
        'workgraph.aiida_int',
        'workgraph.aiida_float',
        'workgraph.aiida_string',
        'workgraph.aiida_bool',
    ]:
        if isinstance(value, TaggedValue):
            value = value.__wrapped__
        if value is not None and isinstance(value, orm.Data):
            return value.value
        else:
            return value
    elif isinstance(value, orm.Data):
        # avoid modifying the original attributes
        content = deepcopy(value.backend_entity.attributes)
        content['node_type'] = value.node_type
        return content


def process_properties(task: Dict) -> Dict:
    """Extract raw values."""
    result = {}
    for name, prop in task.get('properties', {}).items():
        identifier = prop['identifier']
        value = prop.get('value')
        result[name] = {
            'identifier': identifier,
            'value': get_raw_value(identifier, value),
        }
    #
    for name, input in task.get('input_sockets', {}).get('sockets', {}).items():
        if input.get('property'):
            prop = input['property']
            identifier = prop['identifier']
            value = prop.get('value')
            result[name] = {
                'identifier': identifier,
                'value': get_raw_value(identifier, value),
            }

    return result


def workgraph_to_short_json(wgdata: Dict[str, Union[str, List, Dict]]) -> Dict[str, Union[str, Dict]]:
    """Export a workgraph to a rete js editor data."""

    wgdata_short = {
        'name': wgdata['name'],
        'uuid': wgdata.get('uuid', ''),
        'state': wgdata.get('state', ''),
        'nodes': {},
        'links': deepcopy(wgdata.get('links', [])),
    }
    #
    for name, task in wgdata['tasks'].items():
        # Add required inputs to nodes
        inputs = []
        for input in task.get('input_sockets', {}).get('sockets', {}).values():
            metadata = input.get('metadata', {}) or {}
            if metadata.get('required', False):
                inputs.append({'name': input['name'], 'identifier': input['identifier']})

        properties = process_properties(task)
        wgdata_short['nodes'][name] = {
            'identifier': task['identifier'],
            'label': task['name'],
            'node_type': task['spec']['node_type'].upper(),
            'inputs': inputs,
            'properties': properties,
            'outputs': [],
            'position': task.get('position', [0, 0]),
            'children': task.get('children', []),
        }

    # Add links to nodes
    for link in wgdata_short.get('links', []):
        wgdata_short['nodes'][link['to_node']]['inputs'].append(
            {
                'name': link['to_socket'],
            }
        )
        wgdata_short['nodes'][link['from_node']]['outputs'].append(
            {
                'name': link['from_socket'],
            }
        )

    # remove the inputs socket of "graph_inputs"
    if 'graph_inputs' in wgdata_short['nodes']:
        wgdata_short['nodes']['graph_inputs']['inputs'] = []
    # remove the empty graph-level nodes
    for name in ['graph_inputs', 'graph_outputs', 'graph_ctx']:
        if name in wgdata_short['nodes']:
            node = wgdata_short['nodes'][name]
            if len(node['inputs']) == 0 and len(node['outputs']) == 0:
                del wgdata_short['nodes'][name]

    return wgdata_short


def wait_to_link(wgdata: Dict[str, Any]) -> None:
    """Convert wait attribute to link."""
    for name, task in wgdata['tasks'].items():
        for wait_task in task['wait']:
            if wait_task in wgdata['tasks']:
                wgdata['links'].append(
                    {
                        'from_node': wait_task,
                        'from_socket': '_wait',
                        'to_node': name,
                        'to_socket': '_wait',
                    }
                )


def shallow_copy_nested_dict(d):
    """Recursively copies only the dictionary structure but keeps value references."""
    from plumpy.utils import AttributesFrozendict

    if isinstance(d, (dict, AttributesFrozendict)):
        return {key: shallow_copy_nested_dict(value) for key, value in d.items()}
    return d


def make_json_serializable(data):
    """Recursively convert AiiDA objects to JSON-serializable structures."""
    from collections.abc import Mapping, Sequence

    if isinstance(data, orm.Data):
        # Return an int if it's an orm.Int, or a more general dict for other data
        return {
            '__aiida_class__': data.__class__.__name__,
            'uuid': str(data.uuid),
        }
    elif isinstance(data, Mapping):
        return {k: make_json_serializable(v) for k, v in data.items()}
    elif isinstance(data, Sequence) and not isinstance(data, (str, bytes)):
        return [make_json_serializable(item) for item in data]
    else:
        return data


def resolve_tagged_values(inputs: Dict[str, Any]) -> None:
    """Recursively resolve all TaggedValue either in a dictionary or a TaggedValue."""
    if isinstance(inputs, dict):
        for key, value in inputs.items():
            if isinstance(value, TaggedValue):
                inputs[key] = value.__wrapped__
            else:
                resolve_tagged_values(value)


def serialize_graph_level_data(
    input_socket: Dict[str, Any],
    port_schema: SocketSpec | Dict[str, Any],
    serializers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Recursively walk over the sockets and convert raw Python
    values to AiiDA Data nodes, if needed.
    """
    from aiida_pythonjob.utils import serialize_ports

    resolve_tagged_values(input_socket)
    return serialize_ports(
        python_data=input_socket,
        port_schema=port_schema,
        serializers=serializers or {},
    )


def resolve_node_link_managers(data: Any) -> Any:
    """Recursively resolve all NodeLinksManagers either in a dictionary or a NodeLinksManager."""
    if isinstance(data, dict):
        results = {}
        for key, value in data.items():
            results[key] = resolve_node_link_managers(value)
        return results
    elif isinstance(data, orm.NodeLinksManager):
        return convert_node_link_manager_to_dict(data)
    else:
        return data


def convert_node_link_manager_to_dict(
    node_link_manager: orm.NodeLinksManager,
) -> Dict[str, Any]:
    """Recursively convert a NodeLinksManager to a dictionary representation."""
    data = {}
    for name in node_link_manager._get_keys():
        item = node_link_manager._get_node_by_link_label(name)
        if isinstance(item, orm.NodeLinksManager):
            data[name] = convert_node_link_manager_to_dict(item)
        else:
            data[name] = item
    return data


def get_process_summary(node: orm.ProcessNode | int, data: str = ['outputs']) -> None:
    """Get the outputs of a process node."""
    from aiida.common.links import LinkType
    from aiida.cmdline.utils.common import format_nested_links

    node = orm.load_node(node) if isinstance(node, int) else node
    result = ''
    if 'inputs' in data:
        nodes_input = node.base.links.get_incoming(link_type=(LinkType.INPUT_CALC, LinkType.INPUT_WORK))
        result += f'\n{format_nested_links(nodes_input.nested(), headers=["Inputs", "PK", "Type"])}'

    if 'outputs' in data:
        nodes_output = node.base.links.get_outgoing(link_type=(LinkType.CREATE, LinkType.RETURN))
        result += f'\n{format_nested_links(nodes_output.nested(), headers=["Outputs", "PK", "Type"])}'
    return result


def call_depth_from_node(node: str | int | orm.Node) -> int:
    node = orm.load_node(node) if not isinstance(node, orm.Node) else node
    depth = 0
    while getattr(node, 'caller', None) is not None:
        depth += 1
        node = node.caller
    return depth
