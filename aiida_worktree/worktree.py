import scinode
import aiida


class WorkTree(scinode.core.nodetree.NodeTree):
    """
    Represents a working tree for AiiDA's worktree engine. The class extends from scinode's NodeTree and provides
    methods to run, submit tasks, wait for tasks to finish, and update the process status. It is used to handle
    various states of a worktree process and provides convenient operations to interact with it.

    Attributes:
        process (aiida.orm.ProcessNode): The process node that represents the process status and other details.
        state (str): The current state of the worktree process.
        pk (int): The primary key of the process node.
    """

    def __init__(self, name="WorkTree", **kwargs):
        """
        Initialize a WorkTree instance.

        Args:
            name (str, optional): The name of the WorkTree. Defaults to 'WorkTree'.
            **kwargs: Additional keyword arguments to be passed to the NodeTree class.
        """
        super().__init__(name, **kwargs)
        self.ctx = {}
        self.starts = []
        self.is_while = False
        self.conditions = []

    def run(self):
        """
        Run the AiiDA worktree process and update the process status. The method uses AiiDA's engine to run
        the process and then calls the update method to update the state of the process.
        """
        from aiida_worktree.engine.worktree import WorkTree

        ntdata = self.to_dict()
        all = {"nt": ntdata}
        _result, self.process = aiida.engine.run_get_node(WorkTree, **all)
        self.update()

    def submit(self, wait=False, timeout=60):
        """
        Submit the AiiDA worktree process and optionally wait for it to finish.

        Args:
            wait (bool, optional): If True, the function will wait until the process finishes. Defaults to False.
            timeout (int, optional): The maximum time in seconds to wait for the process to finish. Defaults to 60.
        """
        from aiida_worktree.engine.worktree import WorkTree
        from aiida_worktree.utils import merge_properties

        ntdata = self.to_dict()
        merge_properties(ntdata)
        all = {"nt": ntdata}
        self.process = aiida.engine.submit(WorkTree, **all)
        if wait:
            self.wait(timeout=timeout)

    def to_dict(self):
        ntdata = super().to_dict()
        for node in self.nodes:
            ntdata["nodes"][node.name]["to_ctx"] = getattr(node, "to_ctx", [])
        ntdata["ctx"] = self.ctx
        ntdata["starts"] = self.starts
        ntdata["is_while"] = self.is_while
        ntdata["conditions"] = self.conditions
        return ntdata

    def wait(self, timeout=50):
        """
        Periodically checks and waits for the AiiDA worktree process to finish until a given timeout.

        Args:
            timeout (int): The maximum time in seconds to wait for the process to finish. Defaults to 50.
        """
        import time

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
        self.pk = self.process.pk
        outgoing = self.process.base.links.get_outgoing()
        for link in outgoing.all():
            node = link.node
            if isinstance(node, aiida.orm.ProcessNode) and getattr(
                node, "process_state", False
            ):
                self.nodes[link.link_label].state = node.process_state.value.upper()
                self.nodes[link.link_label].node = node
                self.nodes[link.link_label].pk = node.pk
            elif isinstance(node, aiida.orm.Data):
                label = link.link_label.split("__", 1)[1]
                if label in self.nodes.keys():
                    self.nodes[label].state = "FINISHED"
                    self.nodes[label].node = node
                    self.nodes[label].pk = node.pk
