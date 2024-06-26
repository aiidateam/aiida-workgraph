import aiida.orm
import node_graph
import aiida
from aiida.manage import get_manager
from aiida_workgraph.tasks import task_pool
import time
from aiida_workgraph.collection import TaskCollection
from aiida_workgraph.utils.graph import (
    task_deletion_hook,
    task_creation_hook,
    link_creation_hook,
    link_deletion_hook,
)
from aiida_workgraph.widget import NodeGraphWidget
from typing import Any, Dict, List, Optional


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
        self._widget = NodeGraphWidget(parent=self)

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
        # init a process
        runner = get_manager().get_runner()
        process_inited = WorkGraphEngine(runner=runner, inputs=inputs)
        self.process = process_inited.node
        # save workgraph data into process node
        self.save_to_base(inputs["wg"])
        result = aiida.engine.run(process_inited)
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

    def to_dict(self) -> Dict[str, Any]:
        wgdata = super().to_dict()
        self.context["sequence"] = self.sequence
        # only alphanumeric and underscores are allowed
        wgdata["context"] = {
            key.replace(".", "__"): value for key, value in self.context.items()
        }
        wgdata.update(
            {
                "max_iteration": self.max_iteration,
                "execution_count": self.execution_count,
                "workgraph_type": self.workgraph_type,
                "conditions": self.conditions,
                "max_number_jobs": self.max_number_jobs,
            }
        )
        wgdata["tasks"] = wgdata.pop("nodes")

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

        self.state = self.process.process_state.value.upper()
        outgoing = self.process.base.links.get_outgoing()
        for link in outgoing.all():
            node = link.node
            # the link is added in order
            # so the restarted node will be the last one
            # thus the task is correct
            if isinstance(node, aiida.orm.ProcessNode) and getattr(
                node, "process_state", False
            ):
                self.tasks[link.link_label].process = node
                self.tasks[link.link_label].state = node.process_state.value.upper()
                self.tasks[link.link_label].node = node
                self.tasks[link.link_label].pk = node.pk
                self.tasks[link.link_label].ctime = node.ctime
                self.tasks[link.link_label].mtime = node.mtime
                if self.tasks[link.link_label].state == "FINISHED":
                    # update the output sockets
                    i = 0
                    for socket in self.tasks[link.link_label].outputs:
                        if self.tasks[link.link_label].node_type == "graph_builder":
                            if not getattr(node.outputs, "group_outputs", False):
                                continue
                            socket.value = getattr(
                                node.outputs.group_outputs, socket.name, None
                            )
                        else:
                            socket.value = getattr(node.outputs, socket.name, None)
                        i += 1
            elif isinstance(node, aiida.orm.Data):
                if link.link_label.startswith(
                    "group_outputs__"
                ) or link.link_label.startswith("new_data__"):
                    label = link.link_label.split("__", 1)[1]
                    if label in self.tasks.keys():
                        self.tasks[label].state = "FINISHED"
                        self.tasks[label].node = node
                        self.tasks[label].pk = node.pk
                elif link.link_label == "execution_count":
                    self.execution_count = node.value
        # read results from the process outputs
        for task in self.tasks:
            if task.node_type.upper() == "DATA":
                if not getattr(self.process.outputs, "new_data", False):
                    continue
                task.outputs[0].value = getattr(
                    self.process.outputs.new_data, task.name, None
                )
            # for normal tasks, we try to read the results from the extras of the task
            # this is disabled for now
            # if task.node_type.upper() == "NORMAL":
            #     results = self.process.base.extras.get(
            #         f"nodes__results__{task.name}", {}
            #     )
            #     for key, value in results.items():
            #         # if value is an AiiDA data node, we don't need to deserialize it
            #         deserializer = node.outputs[key].get_deserialize()
            #         executor = get_executor(deserializer)[0]
            #         try:
            #             value = executor(bytes(value))
            #         except Exception:
            #             pass
            #         node.outputs[key].value = value
        self._widget.states = {task.name: node.state for node in self.tasks}

    @property
    def pk(self) -> Optional[int]:
        return self.process.pk if self.process else None

    @classmethod
    def from_dict(cls, wgdata: Dict[str, Any]) -> "WorkGraph":
        if "tasks" in wgdata:
            wgdata["nodes"] = wgdata.pop("tasks")
        return super().from_dict(wgdata)

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

    def _repr_mimebundle_(self, *args, **kwargs):
        # if ipywdigets > 8.0.0, use _repr_mimebundle_ instead of _ipython_display_
        self._widget.from_workgraph(self)
        if hasattr(self._widget, "_repr_mimebundle_"):
            return self._widget._repr_mimebundle_(*args, **kwargs)
        else:
            return self._widget._ipython_display_(*args, **kwargs)

    def to_html(self, output: str = None, **kwargs):
        """Write a standalone html file to visualize the workgraph."""
        self._widget.from_workgraph(self)
        return self._widget.to_html(output=output, **kwargs)
