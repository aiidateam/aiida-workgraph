import aiida.orm
import node_graph
import aiida
from aiida_workgraph.nodes import node_pool
import time
from aiida_workgraph.collection import WorkGraphNodeCollection
from aiida_workgraph.utils.graph import (
    node_deletion_hook,
    node_creation_hook,
    link_creation_hook,
    link_deletion_hook,
)
from aiida_workgraph.widget import NodeGraphWidget
from typing import Any, Dict, List, Optional


class WorkGraph(node_graph.NodeGraph):
    """Build a node-based workflow AiiDA's workgraph engine.

    The class extends from NodeGraph and provides methods to run,
    submit tasks, wait for tasks to finish, and update the process status.
    It is used to handle various states of a workgraph process and provides
    convenient operations to interact with it.

    Attributes:
        process (aiida.orm.ProcessNode): The process node that represents the process status and other details.
        state (str): The current state of the workgraph process.
        pk (int): The primary key of the process node.
    """

    node_pool = node_pool

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
        self.nodes = WorkGraphNodeCollection(self, pool=self.node_pool)
        self.nodes.post_deletion_hooks = [node_deletion_hook]
        self.nodes.post_creation_hooks = [node_creation_hook]
        self.links.post_creation_hooks = [link_creation_hook]
        self.links.post_deletion_hooks = [link_deletion_hook]
        self._widget = NodeGraphWidget(parent=self)

    def run(self, inputs: Optional[Dict[str, Any]] = None) -> Any:
        """
        Run the AiiDA workgraph process and update the process status. The method uses AiiDA's engine to run
        the process and then calls the update method to update the state of the process.
        """
        from aiida_workgraph.engine.workgraph import WorkGraph as WorkGraphEngine
        from aiida_workgraph.utils import merge_properties
        from aiida.manage import manager

        # set node inputs
        if inputs is not None:
            for name, input in inputs.items():
                self.nodes[name].set(input)
        # One can not run again if the process is alreay created. otherwise, a new process node will
        # be created again.
        if self.process is not None:
            print("Your workgraph is already created. Please use the submit() method.")
            return
        wgdata = self.to_dict()
        merge_properties(wgdata)
        inputs = {"wg": wgdata}
        # init a process
        runner = manager.get_manager().get_runner()
        process_inited = WorkGraphEngine(runner=runner, inputs=inputs)
        self.process = process_inited.node
        # save workgraph data into process node
        self.save_to_base(wgdata)
        result = aiida.engine.run(process_inited)
        self.update()
        return result

    def submit(
        self,
        inputs: Optional[Dict[str, Any]] = None,
        wait: bool = False,
        timeout: int = 60,
        restart: bool = False,
        new: bool = False,
    ) -> aiida.orm.ProcessNode:
        """Submit the AiiDA workgraph process and optionally wait for it to finish.
        Args:
            wait (bool): Wait for the process to finish.
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 60.
            restart (bool): Restart the process, and reset the modified nodes, then only re-run the modified nodes.
            new (bool): Submit a new process.
        """
        from aiida.manage import get_manager

        # set node inputs
        if inputs is not None:
            for name, input in inputs.items():
                self.nodes[name].set(input)

        process_controller = get_manager().get_process_controller()
        # Create a new submission
        if self.process is not None and new:
            self.reset()
        # Create a restart submission
        # save the current process node as restart_process
        # so that the WorkGraphSaver can compare the difference, and reset the modified nodes
        if restart:
            self.restart_process = self.process
            self.process = None
        # save the workgraph to the process node
        self.save()
        if self.process.process_state.value.upper() not in ["CREATED"]:
            return "Error!!! The process has already been submitted and finished."
        # launch the process, send the task to RabbitMA
        # TODO in case of "[ERROR] Process<3705> is unreachable."
        process_controller.continue_process(self.process.pk)
        if wait:
            self.wait(timeout=timeout)
        return self.process

    def save(self, metadata: Optional[Dict[str, Any]] = None) -> None:
        """Save the udpated workgraph to the process
        This is only used for a running workgraph.
        Save the AiiDA workgraph process and update the process status.
        """
        from aiida_workgraph.engine.workgraph import WorkGraph as WorkGraphEngine
        from aiida_workgraph.utils import merge_properties

        wgdata = self.to_dict()
        merge_properties(wgdata)
        metadata = metadata or {}
        inputs = {"wg": wgdata, "metadata": metadata}
        if self.process is None:
            # init a process node
            process_inited = WorkGraphEngine(inputs=inputs)
            process_inited.runner.persister.save_checkpoint(process_inited)
            self.process = process_inited.node
            self.process_inited = process_inited
            print(f"WorkGraph node created, PK: {self.process.pk}")
        self.save_to_base(wgdata)
        self.update()

    def save_to_base(self, wgdata: Dict[str, Any]) -> None:
        """Save new wgdata to base.extras.
        It will first check the difference, and reset nodes if needed.
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
        of the nodes that are outgoing from the process node. This includes updating the state of process nodes
        linked to the current process, and data nodes linked to the current process.
        """
        # from aiida_workgraph.utils import get_executor

        self.state = self.process.process_state.value.upper()
        outgoing = self.process.base.links.get_outgoing()
        for link in outgoing.all():
            node = link.node
            # the link is added in order
            # so the restarted node will be the last one
            # thus the node is correct
            if isinstance(node, aiida.orm.ProcessNode) and getattr(
                node, "process_state", False
            ):
                self.nodes[link.link_label].process = node
                self.nodes[link.link_label].state = node.process_state.value.upper()
                self.nodes[link.link_label].node = node
                self.nodes[link.link_label].pk = node.pk
                self.nodes[link.link_label].ctime = node.ctime
                self.nodes[link.link_label].mtime = node.mtime
                if self.nodes[link.link_label].state == "FINISHED":
                    # update the output sockets
                    for socket in self.nodes[link.link_label].outputs:
                        if self.nodes[link.link_label].node_type == "graph_builder":
                            if not getattr(node.outputs, "group_outputs", False):
                                continue
                            socket.value = getattr(
                                node.outputs.group_outputs, socket.name, None
                            )
                        else:
                            socket.value = getattr(node.outputs, socket.name, None)
            elif isinstance(node, aiida.orm.Data):
                if link.link_label.startswith(
                    "group_outputs__"
                ) or link.link_label.startswith("new_data__"):
                    label = link.link_label.split("__", 1)[1]
                    if label in self.nodes.keys():
                        self.nodes[label].state = "FINISHED"
                        self.nodes[label].node = node
                        self.nodes[label].pk = node.pk
                elif link.link_label == "execution_count":
                    self.execution_count = node.value
        # read results from the process outputs
        for node in self.nodes:
            if node.node_type.upper() == "DATA":
                if not getattr(self.process.outputs, "new_data", False):
                    continue
                node.outputs[0].value = getattr(
                    self.process.outputs.new_data, node.name, None
                )
            # for normal nodes, we try to read the results from the extras of the node
            # this is disabled for now
            # if node.node_type.upper() == "NORMAL":
            #     results = self.process.base.extras.get(
            #         f"nodes__results__{node.name}", {}
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
        self._widget.states = {node.name: node.state for node in self.nodes}

    @property
    def pk(self) -> Optional[int]:
        return self.process.pk if self.process else None

    @classmethod
    def load(cls, pk: int) -> Optional["WorkGraph"]:
        """
        Load the process node with the given primary key.

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
        for node in self.nodes:
            table.append([node.name, node.pk, node.state])
        print("-" * 80)
        print("WorkGraph: {}, PK: {}, State: {}".format(self.name, self.pk, self.state))
        print("-" * 80)
        # show nodes
        print("Nodes:")
        print(tabulate(table, headers=["Name", "PK", "State"]))
        print("-" * 80)

    def pause(self) -> None:
        """Pause the workgraph."""
        # from aiida.engine.processes import control
        # try:
        # control.pause_processes([self.process])
        import os

        os.system("verdi process pause {}".format(self.process.pk))

    def pause_nodes(self, nodes: List[str]) -> None:
        """
        Pause the given nodes
        """

    def play_nodes(self, nodes: List[str]) -> None:
        """
        Play the given nodes
        """

    def reset(self) -> None:
        """Reset the workgraph."""

        self.process = None
        for node in self.nodes:
            node.reset()
        self.sequence = []
        self.conditions = []
        self.context = {}
        self.state = "CREATED"

    def extend(self, wg: "WorkGraph", prefix: str = "") -> None:
        """Append a workgraph to the current workgraph.
        prefix is used to add a prefix to the node names.
        """
        for node in wg.nodes:
            node.name = prefix + node.name
            node.wait = [prefix + w for w in node.wait] if node.wait else []
            node.parent = self
            self.nodes.append(node)
        # self.sequence.extend([prefix + node for node in wg.sequence])
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
