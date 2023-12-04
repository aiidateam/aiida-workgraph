import node_graph
import aiida
from aiida_worktree.nodes import node_pool
import time


class WorkTree(node_graph.NodeGraph):
    """Build a node-based workflow AiiDA's worktree engine.

    The class extends from NodeGraph and provides methods to run,
    submit tasks, wait for tasks to finish, and update the process status.
    It is used to handle various states of a worktree process and provides
    convenient operations to interact with it.

    Attributes:
        process (aiida.orm.ProcessNode): The process node that represents the process status and other details.
        state (str): The current state of the worktree process.
        pk (int): The primary key of the process node.
    """

    node_pool = node_pool

    def __init__(self, name="WorkTree", **kwargs):
        """
        Initialize a WorkTree instance.

        Args:
            name (str, optional): The name of the WorkTree. Defaults to 'WorkTree'.
            **kwargs: Additional keyword arguments to be passed to the WorkTree class.
        """
        super().__init__(name, **kwargs)
        self.ctx = {}
        self.worktree_type = "NORMAL"
        self.sequence = []
        self.conditions = []
        self.process = None

    def run(self):
        """
        Run the AiiDA worktree process and update the process status. The method uses AiiDA's engine to run
        the process and then calls the update method to update the state of the process.
        """
        from aiida_worktree.engine.worktree import WorkTree as WorkTreeEngine
        from aiida_worktree.utils import merge_properties
        from aiida.orm.utils.serialize import serialize
        from aiida.manage import manager

        # One can not run again if the process is alreay created. otherwise, a new process node will
        # be created again.
        if self.process is not None:
            print("Your worktree is already created. Please use the submit() method.")
            return
        wtdata = self.to_dict()
        merge_properties(wtdata)
        inputs = {"worktree": wtdata}
        # init a process
        runner = manager.get_manager().get_runner()
        process_inited = WorkTreeEngine(runner=runner, inputs=inputs)
        self.process = process_inited.node
        # save worktree data into process node
        self.process.base.extras.set("worktree", serialize(wtdata))
        result = aiida.engine.run(process_inited)
        self.update()
        return result

    def submit(self, wait=False, timeout=60):
        """
        Submit the AiiDA worktree process and optionally wait for it to finish.

        Args:
            wait (bool, optional): If True, the function will wait until the process finishes. Defaults to False.
            timeout (int, optional): The maximum time in seconds to wait for the process to finish. Defaults to 60.
        """
        from aiida.engine.processes import control

        self.save()
        if self.process.process_state.value.upper() not in ["CREATED"]:
            return "Error!!! The process has already been submitted and finished."
        # TODO in case of "[ERROR] Process<3705> is unreachable."
        control.play_processes([self.process])
        if wait:
            self.wait(timeout=timeout)

    def save(self):
        """Save the udpated worktree to the process
        This is only used for a running worktree.
        Save the AiiDA worktree process and update the process status.
        """
        from aiida_worktree.engine.worktree import WorkTree as WorkTreeEngine
        from aiida_worktree.utils import merge_properties
        from aiida.manage import manager

        wtdata = self.to_dict()
        merge_properties(wtdata)
        inputs = {"worktree": wtdata}
        runner = manager.get_manager().get_runner()
        if self.process is None:
            # init a process node
            process_inited = WorkTreeEngine(runner=runner, inputs=inputs)
            runner.persister.save_checkpoint(process_inited)
            # return the future result
            # future = runner.controller.continue_process(process_inited.pid)
            self.process = process_inited.node
            # start = time.time()
            # while not future.done():
            #     time.sleep(1)
            #     if time.time() - start > 5:
            #         print("Worktree dosen't save properly.")
            #         return
            # print(f"WorkTree node crated, PK: {self.process.pk}")
        self.save_to_base(wtdata)
        self.update()

    def save_to_base(self, wtdata):
        """Save new wtdata to base.extras.
        It will first check the difference, and reset nodes if needed.
        """
        from aiida.orm.utils.serialize import serialize

        self.process.base.extras.set("worktree", serialize(wtdata))

    def to_dict(self):
        wtdata = super().to_dict()
        self.ctx["sequence"] = self.sequence
        # only alphanumeric and underscores are allowed
        wtdata["ctx"] = {
            key.replace(".", "__"): value for key, value in self.ctx.items()
        }
        wtdata["worktree_type"] = self.worktree_type
        wtdata["conditions"] = self.conditions

        return wtdata

    def wait(self, timeout=50):
        """
        Periodically checks and waits for the AiiDA worktree process to finish until a given timeout.

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

    def update(self):
        """
        Update the current state and primary key of the process node as well as the state, node and primary key
        of the nodes that are outgoing from the process node. This includes updating the state of process nodes
        linked to the current process, and data nodes linked to the current process.
        """
        self.state = self.process.process_state.value.upper()
        outgoing = self.process.base.links.get_outgoing()
        for link in outgoing.all():
            node = link.node
            if isinstance(node, aiida.orm.ProcessNode) and getattr(
                node, "process_state", False
            ):
                self.nodes[link.link_label].process = node
                self.nodes[link.link_label].state = node.process_state.value.upper()
                self.nodes[link.link_label].node = node
                self.nodes[link.link_label].pk = node.pk
            elif isinstance(node, aiida.orm.Data):
                label = link.link_label.split("__", 1)[1]
                if label in self.nodes.keys():
                    self.nodes[label].state = "FINISHED"
                    self.nodes[label].node = node
                    self.nodes[label].pk = node.pk

    @property
    def pk(self):
        return self.process.pk if self.process else None

    @classmethod
    def load(cls, pk):
        """
        Load the process node with the given primary key.

        Args:
            pk (int): The primary key of the process node.
        """
        from aiida.orm.utils.serialize import deserialize_unsafe

        process = aiida.orm.load_node(pk)
        wtdata = deserialize_unsafe(process.base.extras.get("worktree"))
        wt = cls.from_dict(wtdata)
        wt.process = process
        wt.update()
        return wt

    def show(self):
        """
        Print the current state of the worktree process.
        """
        from tabulate import tabulate

        table = []
        self.update()
        for node in self.nodes:
            table.append([node.name, node.pk, node.state])
        print("-" * 80)
        print("WorkTree: {}, PK: {}, State: {}".format(self.name, self.pk, self.state))
        print("-" * 80)
        # show nodes
        print("Nodes:")
        print(tabulate(table, headers=["Name", "PK", "State"]))
        print("-" * 80)

    def pause_nodes(self, nodes):
        """
        Pause the given nodes
        """

    def play_nodes(self, nodes):
        """
        Play the given nodes
        """
