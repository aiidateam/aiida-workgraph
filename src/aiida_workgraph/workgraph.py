from __future__ import annotations

import aiida.orm
from aiida_workgraph.tasks.builtins import ContextTask
from aiida_workgraph.sockets.builtins import SocketAny
import node_graph
import aiida
import aiida.orm
from node_graph.link import NodeLink
from node_graph.socket import NodeSocketNamespace, NodeSocket

from aiida_workgraph.socket import TaskSocket
from aiida_workgraph.tasks import task_pool
from aiida_workgraph.task import Task
import time
from aiida_workgraph.collection import TaskCollection
from aiida_workgraph.utils.graph import (
    task_deletion_hook,
    task_creation_hook,
    link_creation_hook,
    link_deletion_hook,
)
from aiida_workgraph.orm.mapping import type_mapping
from aiida_workgraph.property import ContextProperty
from typing import Any, Dict, List, Optional, Union

from node_graph_widget import NodeGraphWidget


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

    node_pool = task_pool
    
    class ContexTaskWrapper:
        """Represent a namespace of a Context in the AiiDA WorkGraph."""

        _identifier = "workgraph.namespace"
        _socket_property_class = ContextProperty
        _type_mapping: dict = type_mapping

        def __init__(self, name: str, workgraph: 'WorkGraph', *args, **kwargs):
            # It is a bit tricky because we need to set attributes like the context task
            self._initialization = True
            # I did now ownership to have more flexibility, not sure yet what pattern makes more sense
            #self.nsn = NodeSocketNamespace(name, *args, entry_point="aiida_workgraph.socket", **kwargs) # TODO why do I need  entry_point again?
            self._workgraph = workgraph
            self._context_task = self._workgraph.add_task(ContextTask, name="_CONTEXT")
            self._variable_instances = {}
            self._initialization = False

        #def _new(...): 
        #    # TODO we need at some point replace this to aiida valid labels, not sure where
        #    #wgdata["context"] = {
        #    #    key.replace(".", "__"): value for key, value in self.context.items()
        #    #}

        def __setattr__(self, name: str, value: Any, /) -> None:
            if name == "_initialization" or self._initialization:
                return super().__setattr__(name, value)

            if name.startswith("_"):
                raise ValueError("Socket names starting with an underscores are only for internal usage and are therefore not allowed.")
                
            if name not in self._variable_instances:
                self._variable_instances[name] = 0
            else:
                self._variable_instances[name] += 1

            variable_name = f"{name}{self._variable_instances[name]}"

            # TODO change link limit to some constant constant
            socket = self._context_task.add_variable(variable_name, link_limit=1000)
            #socket = self.nsn._new("workgraph.any", name, link_limit=1000)
            #socket._node = self._context_on 

            if isinstance(value, TaskSocket):
                # purely GUI, the value is populated in engine
                #  socket                 value         
                # wg.context.x = add_task.outputs.result
                if isinstance(value._node, Task):
                    self._workgraph.add_link(value, socket)
                    # TODO(important) to implement this we need multiple tasks otherwise they will all wait for one
                    #self._context_task.waiting_on.add(value._node)
                else:
                    raise TypeError(f"When setting context variable, got node type {type(value._node)} but need Task node.")
            elif isinstance(value, aiida.orm.Data):
                #  socket         value         
                # wg.context.x = orm.Int(5)
                socket.value = value
                # TODO add case where it is convertible to an int type
            else:
                raise TypeError(f"value {value} has invalid type {type(value)}.")

            
        def __getattr__(self, name: str, /) -> Any:
            # Assuming scenario
            # add_task(add, x=wg.context.x)
            # Is this dangerous to set and get different sockets? We have no choice
            # TODO need to validate somehow that socket name does not collide with internal names  
            if name in self._variable_instances:
                variable_name = f"{name}{self._variable_instances[name]}"
                return self._context_task.variables[variable_name]
            # TODO raise better error message? "Socket could not be found"
            return super().__getattribute__(name)

        def update_sockets(self, namespace: 'ContextSocketNamespace'):
            raise NotImplementedError("not properly implemented")
            for socket in self.nsn._sockets:
                if (internal_socket := self._sockets.get(socket._name, None)) is None:
                    self._sockets.append(socket)
                else:
                    internal_socket.value = socket.value
                    internal_socket._node = socket._node

    def __init__(self, name: str = "WorkGraph", **kwargs) -> None:
        """
        Initialize a WorkGraph instance.

        Args:
            name (str, optional): The name of the WorkGraph. Defaults to 'WorkGraph'.
            **kwargs: Additional keyword arguments to be passed to the WorkGraph class.
        """
        super().__init__(name, **kwargs)

        # TODO remove
        self.context = {}
        self.workgraph_type = "NORMAL"
        self.sequence = []
        self.conditions = []
        self.process = None
        self.restart_process = None
        self.max_number_jobs = 1000000
        self.execution_count = 0
        self.max_iteration = 1000000
        self.nodes = TaskCollection(self, pool=self.node_pool)
        self.nodes.post_deletion_hooks = [task_deletion_hook]
        self.nodes.post_creation_hooks = [task_creation_hook]
        self.links.post_creation_hooks = [link_creation_hook]
        self.links.post_deletion_hooks = [link_deletion_hook]
        self._error_handlers = {}
        self._widget = NodeGraphWidget(parent=self)

        self._variables = None

    def __setattr__(self, name: str, value: Any, /) -> None:
        return super().__setattr__(name, value)

    @property
    def variables(self) -> NodeSocketNamespace:
        if self._variables is None:
            self._variables = self.ContexTaskWrapper(name="variables", workgraph=self)
        return self._variables

    @property
    def tasks(self) -> TaskCollection:
        """Add alias to `nodes` for WorkGraph"""
        return self.nodes

    def prepare_inputs(
        self, metadata: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        from aiida_workgraph.utils import serialize_workgraph_inputs

        wgdata = self.to_dict()
        serialize_workgraph_inputs(wgdata)
        metadata = metadata or {}
        inputs = {"wg": wgdata, "metadata": metadata}
        return inputs

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
            for name, input in inputs.items():
                if name not in self.tasks:
                    raise KeyError(f"Task {name} not found in WorkGraph.")
                self.tasks[name].set(input)
        # One can not run again if the process is alreay created. otherwise, a new process node will
        # be created again.
        if self.process is not None:
            print("Your workgraph is already created. Please use the submit() method.")
            return
        inputs = self.prepare_inputs(metadata=metadata)
        result, node = aiida.engine.run_get_node(WorkGraphEngine, inputs=inputs)
        self.process = node
        self.update()
        return result

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
            for name, input in inputs.items():
                if name not in self.tasks:
                    raise KeyError(f"Task {name} not found in WorkGraph.")
                self.tasks[name].set(input)

        # save the workgraph to the process node
        self.save(metadata=metadata)
        if self.process.process_state.value.upper() not in ["CREATED"]:
            raise ValueError(f"Process {self.process.pk} has already been submitted.")
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

        inputs = self.prepare_inputs(metadata)
        if self.process is None:
            runner = manager.get_manager().get_runner()
            # init a process node
            process_inited = instantiate_process(runner, WorkGraphEngine, **inputs)
            process_inited.runner.persister.save_checkpoint(process_inited)
            self.process = process_inited.node
            self.process_inited = process_inited
            process_inited.close()
            print(f"WorkGraph process created, PK: {self.process.pk}")
        else:
            self.save_to_base(inputs["wg"])
        self.update()

    def save_to_base(self, wgdata: Dict[str, Any]) -> None:
        """Save new wgdata to base.extras.
        It will first check the difference, and reset tasks if needed.
        """
        from aiida_workgraph.utils.analysis import WorkGraphSaver

        saver = WorkGraphSaver(
            self.process, wgdata, restart_process=self.restart_process
        )
        saver.save()

    def to_dict(self, store_nodes=False) -> Dict[str, Any]:
        from aiida_workgraph.utils import store_nodes_recursely

        wgdata = super().to_dict()
        # save the sequence and context
        # TODO(session) what is this sequence?
        self.context["_sequence"] = self.sequence
        # only alphanumeric and underscores are allowed
        wgdata["context"] = {
            key.replace(".", "__"): value for key, value in self.context.items()
        }
        wgdata.update(
            {
                "restart_process": aiida.orm.Int(self.restart_process.pk)
                if self.restart_process
                else None,
                "max_iteration": self.max_iteration,
                "execution_count": self.execution_count,
                "workgraph_type": self.workgraph_type,
                "conditions": self.conditions,
                "max_number_jobs": self.max_number_jobs,
            }
        )
        # save error handlers
        wgdata["error_handlers"] = self.get_error_handlers()
        wgdata["tasks"] = wgdata.pop("nodes")
        if store_nodes:
            store_nodes_recursely(wgdata)
        return wgdata

    def get_error_handlers(self) -> Dict[str, Any]:
        """Get the error handlers."""
        from aiida.engine import ExitCode

        from aiida_workgraph.utils import build_callable

        error_handlers = {}
        for name, error_handler in self._error_handlers.items():
            error_handler["handler"] = build_callable(error_handler["handler"])
            error_handlers[name] = error_handler
        # convert exit code label (str) to status (int)
        for handler in error_handlers.values():
            for task in handler["tasks"].values():
                exit_codes = []
                for exit_code in task["exit_codes"]:
                    if isinstance(exit_code, int):
                        exit_codes.append(exit_code)
                    elif isinstance(exit_code, ExitCode):
                        exit_codes.append(exit_code.status)
                    else:
                        raise ValueError(
                            f"Exit code {exit_code} is not a valid exit code."
                        )
                task["exit_codes"] = exit_codes
        return error_handlers

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
            "KILLED",
            "PAUSED",
            "FINISHED",
            "FAILED",
            "CANCELLED",
            "EXCEPTED",
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
                print(f"Process {self.process.pk} finished with state: {self.state}")
                return

            time.sleep(interval)

            if time.time() - start > timeout:
                raise TimeoutError(
                    f"Timeout reached after {timeout} seconds while waiting for the WorkGraph: {self.process.pk}. "
                )

    def update(self) -> None:
        """
        Update the current state and primary key of the process node as well as the state, node and primary key
        of the tasks that are outgoing from the process node. This includes updating the state of process nodes
        linked to the current process, and data nodes linked to the current process.
        """
        # from aiida_workgraph.utils import get_executor
        from aiida_workgraph.utils import get_nested_dict, get_processes_latest

        if self.process is None:
            return

        self.state = self.process.process_state.value.upper()
        processes_data = get_processes_latest(self.pk)
        for name, data in processes_data.items():
            self.tasks[name].state = data["state"]
            self.tasks[name].ctime = data["ctime"]
            self.tasks[name].mtime = data["mtime"]
            self.tasks[name].pk = data["pk"]
            if data["pk"] is not None:
                node = aiida.orm.load_node(data["pk"])
                self.tasks[name].process = self.tasks[name].node = node
                if isinstance(node, aiida.orm.ProcessNode) and getattr(
                    node, "process_state", False
                ):
                    # if the node is finished ok, update the output sockets
                    # note the task.state may not be the same as the node.process_state
                    # for example, task.state can be `SKIPPED` if it is inside a conditional block,
                    # even if the node.is_finished_ok is True
                    if node.is_finished_ok:
                        # update the output sockets
                        for socket in self.tasks[name].outputs:
                            if socket._identifier == "workgraph.namespace":
                                socket._value = get_nested_dict(
                                    node.outputs, socket._name, default=None
                                )
                            else:
                                socket.value = get_nested_dict(
                                    node.outputs, socket._name, default=None
                                )
                # read results from the process outputs
                elif isinstance(node, aiida.orm.Data):
                    self.tasks[name].outputs[0].value = node
        execution_count = getattr(self.process.outputs, "execution_count", None)
        self.execution_count = execution_count if execution_count else 0
        if self._widget is not None:
            states = {name: data["state"] for name, data in processes_data.items()}
            self._widget.states = states

    @property
    def pk(self) -> Optional[int]:
        return self.process.pk if self.process else None

    @classmethod
    def from_dict(cls, wgdata: Dict[str, Any]) -> "WorkGraph":

        if "tasks" in wgdata:
            wgdata["nodes"] = wgdata.pop("tasks")
        wg = super().from_dict(wgdata)
        for key in [
            "max_iteration",
            "execution_count",
            "workgraph_type",
            "conditions",
            "max_number_jobs",
        ]:
            if key in wgdata:
                setattr(wg, key, wgdata[key])
        if "error_handlers" in wgdata:
            wg._error_handlers = wgdata["error_handlers"]
        return wg

    @classmethod
    def from_yaml(cls, filename: str = None, string: str = None) -> "WorkGraph":
        """Build WrokGraph from yaml file."""
        import yaml
        from node_graph.utils import yaml_to_dict

        # load data
        if filename:
            with open(filename, "r") as f:
                wgdata = yaml.safe_load(f)
        elif string:
            wgdata = yaml.safe_load(string)
        else:
            raise Exception("Please specific a filename or yaml string.")
        wgdata["nodes"] = wgdata.pop("tasks")
        wgdata = yaml_to_dict(wgdata)
        nt = cls.from_dict(wgdata)
        return nt

    @classmethod
    def load(cls, pk: int) -> Optional["WorkGraph"]:
        """
        Load WorkGraph from the process node with the given primary key.

        Args:
            pk (int): The primary key of the process node.
        """
        from aiida_workgraph.utils import get_workgraph_data

        process = aiida.orm.load_node(pk)
        wgdata = get_workgraph_data(process)
        if wgdata is None:
            raise ValueError(f"WorkGraph data not found for process {pk}.")
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
        print("-" * 80)
        print("WorkGraph: {}, PK: {}, State: {}".format(self.name, self.pk, self.state))
        print("-" * 80)
        print("Tasks:")
        print(tabulate(table, headers=["Name", "PK", "State"]))
        print("-" * 80)

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
                self.tasks[name].action = "PAUSE"
        else:
            _, msg = pause_tasks(self.process.pk, tasks)

        return "Send message to pause tasks."

    def play_tasks(self, tasks: List[str]) -> None:
        """Play the given tasks"""

        from aiida_workgraph.utils.control import play_tasks

        if self.process is None:
            for name in tasks:
                self.tasks[name].action = ""
        else:
            _, msg = play_tasks(self.process.pk, tasks)
        return "Send message to play tasks."

    def kill_tasks(self, tasks: List[str]) -> None:
        """Kill the given tasks"""

        from aiida_workgraph.utils.control import kill_tasks

        if self.process is None:
            for name in tasks:
                self.tasks[name].action = "KILL"
        else:
            _, msg = kill_tasks(self.process.pk, tasks)
        return "Send message to kill tasks."

    def continue_process(self):
        """Continue a saved process by sending the task to RabbitMA.
        Use with caution, this may launch duplicate processes."""
        from aiida.manage import get_manager

        process_controller = get_manager().get_process_controller()
        process_controller.continue_process(self.pk)

    def play(self):
        import os

        os.system("verdi process play {}".format(self.process.pk))

    def restart(self):
        """Create a restart submission."""
        if self.process is None:
            raise ValueError(
                "No process found. One can not restart from a non-existing process."
            )
        # save the current process node as restart_process
        # so that the WorkGraphSaver can compare the difference, and reset the modified tasks
        self.restart_process = self.process
        self.process = None

    def reset(self) -> None:
        """Reset the workgraph to create a new submission."""

        self.process = None
        for task in self.tasks:
            task.reset()
        self.sequence = []
        self.conditions = []
        self.context = {}
        self._variables = None 
        self.state = "PLANNED"

    def extend(self, wg: "WorkGraph", prefix: str = "") -> None:
        """Append a workgraph to the current workgraph.
        prefix is used to add a prefix to the task names.
        """
        for task in wg.tasks:
            task.name = prefix + task.name
            task.parent = self
            self.tasks._append(task)
        # self.sequence.extend([prefix + task for task in wg.sequence])
        # self.conditions.extend(wg.conditions)
        self.context.update(wg.context)
        # TODO _sockets and update sockets need to be public API
        self.variables.update_sockets(wg.variables)
        # links
        for link in wg.links:
            self.links._append(link)

    @property
    def error_handlers(self) -> Dict[str, Any]:
        """Get the error handlers."""
        return self._error_handlers

    def add_error_handler(self, handler, name, tasks: dict = None) -> None:
        """Attach an error handler to the workgraph."""
        self._error_handlers[name] = {"handler": handler, "tasks": tasks}

    def add_task(
        self, identifier: Union[str, callable], name: str = None, **kwargs
    ) -> Task:
        """Add a task to the workgraph."""
        node = self.tasks._new(identifier, name, **kwargs)
        return node

    def add_link(self, source: TaskSocket | Task, target: TaskSocket) -> NodeLink:
        """Add a link between two tasks."""
        if isinstance(source, Task):
            source = source.outputs["_outputs"]
        link = self.links._new(source, target)
        return link

    def to_widget_value(self) -> Dict[str, Any]:
        """Convert the workgraph to a dictionary that can be used by the widget."""
        from aiida_workgraph.utils import workgraph_to_short_json, wait_to_link

        wgdata = self.to_dict()
        wait_to_link(wgdata)
        wgdata = workgraph_to_short_json(wgdata)
        return wgdata

    def _repr_mimebundle_(self, *args, **kwargs):
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.value = self.to_widget_value()
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)

    def to_html(self, output: str = None, **kwargs):
        """Write a standalone html file to visualize the workgraph."""
        self._widget.value = self.to_widget_value()
        return self._widget.to_html(output=output, **kwargs)
