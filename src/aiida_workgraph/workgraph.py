from __future__ import annotations

import node_graph
import aiida
from aiida_workgraph.task import Task
import time
from typing import Any, Dict, List, Optional, Union
from .registry import RegistryHub, registry_hub
from node_graph.analysis import NodeGraphAnalysis
from node_graph.config import BUILTIN_NODES
from node_graph.collection import NodeCollection
from node_graph.socket import BaseSocket, NodeSocketNamespace
from aiida_workgraph.socket_spec import SocketSpecAPI
from node_graph.error_handler import ErrorHandlerSpec


class WorkGraph(node_graph.NodeGraph):
    """Build flexible workflows with AiiDA.

    The class extends from NodeGraph and provides methods to run,
    submit tasks, wait for tasks to finish, and update the process status.
    It is used to handle various states of a workgraph process and provides
    convenient operations to interact with it.

    Attributes:
        process (aiida.orm.ProcessNode): The process node that represents the process status and other details.
        state (str): The current state of the workgraph process.
        pk (int): The primary key of the process node.
    """

    _REGISTRY: Optional[RegistryHub] = registry_hub
    _SOCKET_SPEC_API = SocketSpecAPI

    platform: str = 'aiida_workgraph'

    def __init__(
        self,
        name: str = 'WorkGraph',
        inputs: Optional[type | List[str]] = None,
        outputs: Optional[type | List[str]] = None,
        error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a WorkGraph instance.

        Args:
            name (str, optional): The name of the WorkGraph. Defaults to 'WorkGraph'.
            **kwargs: Additional keyword arguments to be passed to the WorkGraph class.
        """
        super().__init__(name, inputs=inputs, outputs=outputs, **kwargs)
        self.process = None
        self.restart_process = None
        self.max_number_jobs = 1000000
        self.max_iteration = 1000000
        self._error_handlers = error_handlers or {}
        self.analyzer = NodeGraphAnalysis(self)

    @property
    def tasks(self) -> NodeCollection:
        """Add alias to `nodes` for WorkGraph"""
        return self.nodes

    def to_engine_inputs(self, metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        wgdata = self.to_dict(should_serialize=True)
        metadata = metadata or {}
        task_inputs = self.gather_task_inputs(wgdata['tasks'])
        graph_inputs = task_inputs.pop('graph_inputs', {})
        inputs = {
            'metadata': metadata,
            'workgraph_data': wgdata,
            'tasks': task_inputs,
            'graph_inputs': graph_inputs,
        }
        return inputs

    def gather_task_inputs(self, data: Dict[str, Any] = None) -> Dict[str, Any]:
        """Gather the inputs of all tasks."""
        inputs = {}
        for name, task in data.items():
            inputs[name] = task.pop('inputs', {})
        return inputs

    def check_before_run(self) -> bool:
        self.check_required_inputs()
        self.check_modified_tasks()

    def check_required_inputs(self) -> None:
        """Check if all required inputs are provided."""
        missing_inputs = self.find_missing_inputs(self.inputs)
        for task in self.tasks:
            if task.name in BUILTIN_NODES:
                continue
            missing_inputs.extend(self.find_missing_inputs(task.inputs))
        if missing_inputs:
            bullets = '\n'.join(f'  • {p}' for p in sorted(missing_inputs))
            raise ValueError(
                'Missing required inputs:\n'
                f'{bullets}\n\n'
                'How to fix:\n'
                '  1) Provide these values (at build time or by linking from upstream task outputs).\n'
                '  2) If some are intentionally unused, exclude them from the namespace at the call site, e.g.:\n'
                '     Annotated[dict, some_task.inputs, SocketSpecSelect(exclude=["pw.structure", ...])]\n\n'
                "Note: exclude paths are relative to the task's input namespace (e.g. 'pw.structure')."
            )

    def find_missing_inputs(self, socket: BaseSocket) -> List[str]:
        """Check if all required inputs are provided."""
        missing_inputs = []
        for sub_socket in socket:
            if isinstance(sub_socket, NodeSocketNamespace):
                missing_inputs.extend(self.find_missing_inputs(sub_socket))
            else:
                if sub_socket._metadata.required and sub_socket.value is None and len(sub_socket._links) == 0:
                    missing_inputs.append(f'{sub_socket._node.name}.{sub_socket._scoped_name}')
        return missing_inputs

    def check_modified_tasks(self) -> None:
        """Check if there are modified tasks compared to the existing process.
        If there are modified tasks, reset them and their descendants to be re-run.
        """
        existing_process = self._load_existing_process()
        if existing_process:
            diffs = self.analyzer.compare_graphs(existing_process, self)
            self.reset_tasks(diffs['modified_nodes'])

    def run(
        self,
        inputs: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Any:
        """
        Run the AiiDA workgraph process and update the process status. The method uses AiiDA's engine to run
        the process, when the process is finished, update the status of the tasks
        """
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        # set task inputs
        if inputs is not None:
            self.set_inputs(inputs)

        # One can not run again if the process is alreay created. otherwise, a new process node will
        # be created again.
        if self.process is not None:
            raise ValueError(f'Process {self.process.pk} has already been created. Please use the submit() method.')
        self.check_before_run()
        inputs = self.to_engine_inputs(metadata=metadata)
        _, node = aiida.engine.run_get_node(WorkGraphEngine, inputs=inputs)
        self.process = node
        self.update()
        return self.outputs._value

    def submit(
        self,
        inputs: Optional[Dict[str, Any]] = None,
        wait: bool = False,
        timeout: int = 600,
        interval: int = 5,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> aiida.orm.ProcessNode:
        """Submit the AiiDA workgraph process and optionally wait for it to finish.
        Args:
            wait (bool): Wait for the process to finish.
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 600.
            restart (bool): Restart the process, and reset the modified tasks, then only re-run the modified tasks.
            new (bool): Submit a new process.
        """

        # set task inputs
        if inputs is not None:
            self.set_inputs(inputs)

        # save the workgraph to the process node
        self.save(metadata=metadata)
        if self.process.process_state.value.upper() not in ['CREATED']:
            raise ValueError(f'Process {self.process.pk} has already been submitted.')
        self.continue_process()
        # as long as we submit the process, it is a new submission, we should set restart_process to None
        self.restart_process = None
        if wait:
            self.wait(timeout=timeout, interval=interval)
        return self.process

    def save(self, metadata: Optional[Dict[str, Any]] = None) -> None:
        """Save the udpated workgraph to the process
        This is only used for a running workgraph.
        Save the AiiDA workgraph process and update the process status.
        """
        from aiida.manage import manager
        from aiida.engine.utils import instantiate_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        self.check_before_run()
        inputs = self.to_engine_inputs(metadata)
        if self.process is None:
            runner = manager.get_manager().get_runner()
            # init a process node
            process_inited = instantiate_process(runner, WorkGraphEngine, **inputs)
            process_inited.runner.persister.save_checkpoint(process_inited)
            self.process = process_inited.node
            self.process_inited = process_inited
            process_inited.close()
            print(f'WorkGraph process created, PK: {self.process.pk}')
        else:
            self.save_to_base(inputs)
        self.update()

    def save_to_base(self, wgdata: Dict[str, Any]) -> None:
        """Save new wgdata to attribute.
        It will first check the difference, and reset tasks if needed.
        """
        from aiida_workgraph.utils import save_workgraph_data

        save_workgraph_data(self.process, wgdata)

    def _load_existing_process(self):
        """Load an existing workgraph process if available."""
        if self.process:
            return WorkGraph.load(self.process)
        if self.restart_process:
            return WorkGraph.load(self.restart_process)
        return None

    def build_connectivity(self) -> None:
        """Analyze the connectivity of workgraph and save it into dict."""
        connectivity = self.analyzer.build_connectivity()
        return connectivity

    def to_dict(self, include_sockets: bool = False, should_serialize: bool = False) -> Dict[str, Any]:
        """Convert the workgraph to a dictionary."""
        from aiida.orm.utils.serialize import serialize
        from aiida_workgraph.utils import serialize_graph_level_data
        from aiida_pythonjob.data.serializer import all_serializers

        wgdata = super().to_dict(include_sockets=include_sockets, should_serialize=should_serialize)
        wgdata.update(
            {
                'restart_process': self.restart_process.pk if self.restart_process else None,
                'max_iteration': self.max_iteration,
                'max_number_jobs': self.max_number_jobs,
            }
        )
        # save error handlers
        wgdata['error_handlers'] = {name: eh.to_dict() for name, eh in self.get_error_handlers().items()}
        wgdata['tasks'] = wgdata.pop('nodes')
        wgdata['connectivity'] = self.build_connectivity()
        wgdata['process'] = serialize(self.process) if self.process else serialize(None)
        wgdata['metadata']['pk'] = self.process.pk if self.process else None
        if should_serialize:
            # serialize the graph-level tasks
            wgdata['tasks']['graph_inputs']['inputs'] = serialize_graph_level_data(
                wgdata['tasks']['graph_inputs']['inputs'], self.spec.inputs, all_serializers
            )

        return wgdata

    def wait(self, timeout: int = 600, tasks: dict = None, interval: int = 5) -> None:
        """
        Periodically checks and waits for the AiiDA workgraph process to finish until a given timeout.

        Args:
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 600.
            tasks (dict): Optional; specifies task states to wait for in the format {task_name: [acceptable_states]}.
            interval (int): The time interval in seconds between checks. Defaults to 5.

        Raises:
            TimeoutError: If the process does not finish within the given timeout.
        """
        terminating_states = (
            'KILLED',
            'PAUSED',
            'FINISHED',
            'FAILED',
            'CANCELLED',
            'EXCEPTED',
        )
        start = time.time()
        self.update()
        finished = False

        while not finished:
            self.update()

            if tasks is not None:
                states = []
                for name, value in tasks.items():
                    flag = self.tasks[name].state in value
                    states.append(flag)
                finished = all(states)
            else:
                finished = self.state in terminating_states

            if finished:
                print(f'Process {self.process.pk} finished with state: {self.state}')
                return

            time.sleep(interval)

            if time.time() - start > timeout:
                raise TimeoutError(
                    f'Timeout reached after {timeout} seconds while waiting for the WorkGraph: {self.process.pk}. '
                )

    def update(self) -> None:
        """
        Update the current state and primary key of the process node as well as the state, node and primary key
        of the tasks that are outgoing from the process node. This includes updating the state of process nodes
        linked to the current process, and data nodes linked to the current process.
        """
        from aiida_workgraph.utils import get_processes_latest, resolve_node_link_managers

        if self.process is None:
            return

        self.state = self.process.process_state.value.upper()
        processes_data = get_processes_latest(self.pk)
        for name, data in processes_data.items():
            # the mapped tasks are not in the workgraph
            if name not in self.tasks:
                continue
            self.tasks[name].update_state(data)

        if self.widget is not None:
            states = {name: data['state'] for name, data in processes_data.items()}
            self.widget.states = states

        if self.process.is_finished_ok:
            self.outputs._set_socket_value(resolve_node_link_managers(self.process.outputs))

    @property
    def pk(self) -> Optional[int]:
        return self.process.pk if self.process else None

    @classmethod
    def from_dict(cls, wgdata: Dict[str, Any]) -> 'WorkGraph':
        if 'tasks' in wgdata:
            wgdata['nodes'] = wgdata.pop('tasks')
        wg = super().from_dict(wgdata)
        for key in [
            'max_iteration',
            'max_number_jobs',
            'connectivity',
        ]:
            if key in wgdata:
                setattr(wg, key, wgdata[key])
        if 'error_handlers' in wgdata:
            wg._error_handlers = {
                name: ErrorHandlerSpec.from_dict(eh) for name, eh in wgdata.get('error_handlers', {}).items()
            }
        # for zone tasks, add their children
        for task in wg.tasks:
            if hasattr(task, 'children'):
                task.children.add(wgdata['nodes'][task.name].get('children', []))
        # reinstate the tasks
        wgdata['tasks'] = wgdata.pop('nodes')

        return wg

    @classmethod
    def from_yaml(cls, filename: str = None, string: str = None) -> 'WorkGraph':
        """Build WrokGraph from yaml file."""
        import yaml

        # import json
        # from aiida_workgraph.utils import make_json_serializable
        from node_graph.utils import yaml_to_dict

        # import importlib.resources
        # import jsonschema

        if filename:
            with open(filename, 'r') as f:
                wgdata = yaml.safe_load(f)
        elif string:
            wgdata = yaml.safe_load(string)
        else:
            raise Exception('Please specific a filename or yaml string.')
        wgdata['nodes'] = wgdata.pop('tasks')
        wgdata = yaml_to_dict(wgdata)
        wgdata['tasks'] = wgdata.pop('nodes')
        # serialized_data = make_json_serializable(wgdata)
        # with importlib.resources.open_text(
        #     "aiida_workgraph.schemas", "aiida_workgraph.schema.json"
        # ) as f:
        #     schema = json.load(f)
        #     jsonschema.validate(instance=serialized_data, schema=schema)

        nt = cls.from_dict(wgdata)
        return nt

    @classmethod
    def load(cls, pk: int | str | aiida.orm.ProcessNode) -> Optional['WorkGraph']:
        """
        Load WorkGraph from the process node with the given primary key.

        Args:
            pk (int, str, orm.ProcessNode): The primary key or uuid of the process node,
                or the process node itself.
        """
        from aiida_workgraph.orm.workgraph import WorkGraphNode
        from aiida_workgraph.utils import load_workgraph_data

        if isinstance(pk, (int, str)):
            process = aiida.orm.load_node(pk)
        elif isinstance(pk, aiida.orm.ProcessNode):
            process = pk
        else:
            raise ValueError(f'Invalid pk type: {type(pk)}, requires int, str or ProcessNode.')
        if not isinstance(process, WorkGraphNode):
            raise ValueError(f'Process {pk} is not a WorkGraph')
        wgdata = load_workgraph_data(process)
        wg = cls.from_dict(wgdata)
        wg.process = process
        wg.update()
        return wg

    def show(self) -> None:
        """
        Print the current state of the workgraph process.
        """
        from tabulate import tabulate

        table = []
        self.update()
        for task in self.tasks:
            table.append([task.name, task.pk, task.state])
        print('-' * 80)
        print('WorkGraph: {}, PK: {}, State: {}'.format(self.name, self.pk, self.state))
        print('-' * 80)
        print('Tasks:')
        print(tabulate(table, headers=['Name', 'PK', 'State']))
        print('-' * 80)

    # def pause(self) -> None:
    #     """Pause the workgraph."""
    #     from aiida.engine.processes import control
    #     try:
    #         control.pause_processes([self.process])
    #     except Exception as e:
    #         print(f"Pause process failed: {e}")

    def pause_tasks(self, tasks: List[str]) -> None:
        """Pause the given tasks."""
        from aiida_workgraph.utils.control import pause_tasks

        if self.process is None:
            for name in tasks:
                self.tasks[name].action = 'PAUSE'
        else:
            _, msg = pause_tasks(self.process.pk, tasks)

        return 'Send message to pause tasks.'

    def play_tasks(self, tasks: List[str]) -> None:
        """Play the given tasks"""

        from aiida_workgraph.utils.control import play_tasks

        if self.process is None:
            for name in tasks:
                self.tasks[name].action = ''
        else:
            _, msg = play_tasks(self.process.pk, tasks)
        return 'Send message to play tasks.'

    def kill_tasks(self, tasks: List[str]) -> None:
        """Kill the given tasks"""

        from aiida_workgraph.utils.control import kill_tasks

        if self.process is None:
            for name in tasks:
                self.tasks[name].action = 'KILL'
        else:
            _, msg = kill_tasks(self.process.pk, tasks)
        return 'Send message to kill tasks.'

    def reset_tasks(self, tasks: List[str]) -> None:
        from aiida_workgraph.utils.control import reset_tasks

        print(f'Reset tasks: {tasks}')

        if self.process is None:
            for name in tasks:
                self.tasks[name].state = 'PLANNED'
                self.tasks[name].process = None
                child_tasks = self.analyzer.get_all_descendants(self.tasks[name])
                for name in child_tasks:
                    self.tasks[name].state = 'PLANNED'
                    self.tasks[name].process = None
        else:
            _, msg = reset_tasks(self.process.pk, tasks)
        return 'Send message to reset tasks.'

    def continue_process(self):
        """Continue a saved process by sending the task to RabbitMA.
        Use with caution, this may launch duplicate processes."""
        from aiida.manage import get_manager

        process_controller = get_manager().get_process_controller()
        process_controller.continue_process(self.pk)

    def play(self):
        import os

        os.system('verdi process play {}'.format(self.process.pk))

    def restart(self):
        """Create a restart submission."""
        if self.process is None:
            raise ValueError('No process found. One can not restart from a non-existing process.')
        # save the current process node as restart_process
        # so that the WorkGraphSaver can compare the difference, and reset the modified tasks
        self.restart_process = self.process
        self.process = None
        self.state = 'PLANNED'

    def reset(self) -> None:
        """Reset the workgraph to create a new submission."""

        self.process = None
        for task in self.tasks:
            task.reset()
        self.state = 'PLANNED'

    def extend(self, wg: 'WorkGraph', prefix: str = '') -> None:
        """Append a workgraph to the current workgraph.
        prefix is used to add a prefix to the task names.
        """
        for task in wg.tasks:
            # skip the built-in tasks
            # need to fix this in the future
            if task.name in BUILTIN_NODES:
                continue
            task.name = prefix + task.name
            task.graph = self
            self.tasks._append(task)
        self.update_ctx(wg.ctx._value)
        # links
        for link in wg.links:
            # skip the links that are from or to built-in tasks
            if link.from_node.name in BUILTIN_NODES or link.to_node.name in BUILTIN_NODES:
                link.unmount()
                continue
            self.links._append(link)

    def get_error_handlers(self) -> Dict[str, ErrorHandlerSpec]:
        """Get the error handlers."""
        return self._error_handlers

    def add_task(
        self,
        identifier: Union[str, callable],
        name: str = None,
        include_builtins: bool = False,
        **kwargs,
    ) -> Task:
        """Add a task to the workgraph."""
        from aiida_workgraph.decorator import build_task_from_callable
        from aiida.engine import ProcessBuilder
        from aiida_workgraph.utils import get_dict_from_builder
        from aiida_workgraph.tasks.shelljob_task import (
            shelljob,
            _build_shelljob_nodespec,
        )
        from aiida_workgraph.task import Task, TaskHandle
        from aiida_workgraph.tasks.subgraph_task import _build_subgraph_task_nodespec
        from node_graph.node_spec import NodeSpec

        if name in BUILTIN_NODES and not include_builtins:
            raise ValueError(f'Task name {name} can not be used, it is reserved.')

        if isinstance(identifier, str):
            identifier = self._REGISTRY.node_pool[identifier.lower()].load()
        if isinstance(identifier, WorkGraph):
            identifier = _build_subgraph_task_nodespec(identifier, name=name)
        elif isinstance(identifier, ProcessBuilder):
            kwargs = {**kwargs, **get_dict_from_builder(identifier)}
            identifier = build_task_from_callable(identifier.process_class)
        # todo
        elif identifier is shelljob:
            spec = _build_shelljob_nodespec(
                outputs=kwargs.get('outputs'),
                parser_outputs=kwargs.pop('parser_outputs', None),
            )
            identifier = TaskHandle(spec)
        # build the task on the fly if the identifier is a callable
        elif callable(identifier) and not isinstance(identifier, (NodeSpec, TaskHandle, Task)):
            identifier = build_task_from_callable(identifier)
        node = self.tasks._new(identifier, name, **kwargs)
        self._version += 1
        return node

    def to_widget_value(self) -> Dict[str, Any]:
        """Convert the workgraph to a dictionary that can be used by the widget."""
        from aiida_workgraph.utils import workgraph_to_short_json, wait_to_link

        wgdata = self.to_dict(include_sockets=True)
        wait_to_link(wgdata)
        wgdata = workgraph_to_short_json(wgdata)
        return wgdata

    def generate_provenance_graph(self):
        """Generate the provenance graph of the workgraph process."""
        from aiida_workgraph.utils import generate_provenance_graph

        if self.process is None:
            raise ValueError('No process found. Please run or submit the workgraph first.')
        return generate_provenance_graph(self.process.pk)

    def __repr__(self) -> str:
        return f'WorkGraph(name="{self.name}", uuid="{self.uuid}")'

    def __str__(self) -> str:
        return f'WorkGraph(name="{self.name}", uuid="{self.uuid}")'

    def __enter__(self):
        """Called when entering the `with NodeGraph() as ng:` block."""
        from aiida_workgraph.manager import get_current_graph, set_current_graph

        self._previous_graph = get_current_graph()
        set_current_graph(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Called upon leaving the `with NodeGraph() as ng:` block."""
        from aiida_workgraph.manager import set_current_graph

        set_current_graph(self._previous_graph)
        self._previous_graph = None
        return None

    def __call__(self, inputs: Dict[str, Any] = None) -> Any:
        """Call the graph with inputs and return as a task.

        Used in context managers as a simple assignment.

        >>> wg1 = WorkGraph()
        >>> with WorkGraph() as wg2:
        >>>     task_outputs = wg1({'input1': 42, 'input2': 'hello'})
        """
        from aiida_workgraph.manager import get_current_graph

        graph = get_current_graph()
        task = graph.add_task(self)
        inputs = inputs or {}
        task.set_inputs(inputs)
        return task.outputs
