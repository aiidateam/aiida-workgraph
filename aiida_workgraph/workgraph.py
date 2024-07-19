import aiida.orm
import node_graph
import aiida
import node_graph.link
from aiida_workgraph.socket import NodeSocket
from aiida_workgraph import USE_WIDGET
from aiida_workgraph.utils.message import WIDGET_INSTALLATION_MESSAGE
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
from typing import Any, Dict, List, Optional, Union

if USE_WIDGET:
    from aiida_workgraph.widget import NodeGraphWidget


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

    def __init__(self, name: str = "WorkGraph", **kwargs) -> None:
        """
        Initialize a WorkGraph instance.

        Args:
            name (str, optional): The name of the WorkGraph. Defaults to 'WorkGraph'.
            **kwargs: Additional keyword arguments to be passed to the WorkGraph class.
        """
        super().__init__(name, **kwargs)
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
        self.error_handlers = {}
        self._widget = NodeGraphWidget(parent=self) if USE_WIDGET else None

    @property
    def tasks(self) -> TaskCollection:
        """Add alias to `nodes` for WorkGraph"""
        return self.nodes

    def prepare_inputs(self, metadata: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        from aiida_workgraph.utils import (
            merge_properties,
            serialize_pythonjob_properties,
        )

        wgdata = self.to_dict()
        merge_properties(wgdata)
        serialize_pythonjob_properties(wgdata)
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
                if name not in self.tasks.keys():
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
        timeout: int = 60,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> aiida.orm.ProcessNode:
        """Submit the AiiDA workgraph process and optionally wait for it to finish.
        Args:
            wait (bool): Wait for the process to finish.
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 60.
            restart (bool): Restart the process, and reset the modified tasks, then only re-run the modified tasks.
            new (bool): Submit a new process.
        """
        # set task inputs
        if inputs is not None:
            for name, input in inputs.items():
                if name not in self.tasks.keys():
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
            self.wait(timeout=timeout)
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
        import cloudpickle as pickle

        wgdata = super().to_dict()
        self.context["sequence"] = self.sequence
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
        wgdata["error_handlers"] = pickle.dumps(self.error_handlers)
        wgdata["tasks"] = wgdata.pop("nodes")
        if store_nodes:
            for task in wgdata["tasks"].values():
                for prop in task["properties"].values():
                    if isinstance(prop["value"], aiida.orm.Node):
                        prop["value"].store()

        return wgdata

    def wait(self, timeout: int = 50) -> None:
        """
        Periodically checks and waits for the AiiDA workgraph process to finish until a given timeout.

        Args:
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 50.
        """

        start = time.time()
        self.update()
        while self.state not in (
            "KILLED",
            "PAUSED",
            "FINISHED",
            "FAILED",
            "CANCELLED",
            "EXCEPTED",
        ):
            time.sleep(0.5)
            self.update()
            if time.time() - start > timeout:
                return

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
                    if self.tasks[name].state == "FINISHED":
                        # update the output sockets
                        i = 0
                        for socket in self.tasks[name].outputs:
                            socket.value = get_nested_dict(
                                node.outputs, socket.name, allow_none=True
                            )
                            i += 1
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
        import cloudpickle as pickle

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
            wg.error_handlers = pickle.loads(wgdata["error_handlers"])
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
            print("No workgraph data found in the process node.")
            return
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

    def pause(self) -> None:
        """Pause the workgraph."""
        # from aiida.engine.processes import control
        # try:
        # control.pause_processes([self.process])
        import os

        os.system("verdi process pause {}".format(self.process.pk))

    def pause_tasks(self, tasks: List[str]) -> None:
        """Pause the given tasks."""
        from aiida_workgraph.utils.control import pause_tasks

        self.update()

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
        self.state = "PLANNED"

    def extend(self, wg: "WorkGraph", prefix: str = "") -> None:
        """Append a workgraph to the current workgraph.
        prefix is used to add a prefix to the task names.
        """
        for task in wg.tasks:
            task.name = prefix + task.name
            task.wait = [prefix + w for w in task.wait] if task.wait else []
            task.parent = self
            self.tasks.append(task)
        # self.sequence.extend([prefix + task for task in wg.sequence])
        # self.conditions.extend(wg.conditions)
        self.context.update(wg.context)
        # links
        for link in wg.links:
            self.links.append(link)

    def attach_error_handler(self, handler, name, tasks: dict = None) -> None:
        """Attach an error handler to the workgraph."""
        self.error_handlers[name] = {"handler": handler, "tasks": tasks}

    def _repr_mimebundle_(self, *args, **kwargs):

        if self._widget is None:
            print(WIDGET_INSTALLATION_MESSAGE)
            return

        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.from_workgraph(self)
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)

    def add_task(
        self, identifier: Union[str, callable], name: str = None, **kwargs
    ) -> Task:
        """Add a task to the workgraph."""
        node = self.tasks.new(identifier, name, **kwargs)
        return node

    def add_link(
        self, source: NodeSocket, target: NodeSocket
    ) -> node_graph.link.NodeLink:
        """Add a link between two nodes."""
        link = self.links.new(source, target)
        return link

    def to_html(self, output: str = None, **kwargs):
        """Write a standalone html file to visualize the workgraph."""
        if self._widget is None:
            print(WIDGET_INSTALLATION_MESSAGE)
            return
        self._widget.from_workgraph(self)
        return self._widget.to_html(output=output, **kwargs)
