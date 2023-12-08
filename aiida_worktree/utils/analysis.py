import datetime


class WorkTreeSaver:
    """Save a worktree to the database."""

    def __init__(self, process, wtdata, restart_process=None) -> None:
        """Init WorkTreeSaver.

        Args:
            wtdata (dict): data of nodetree to be launched.
        """
        self.process = process
        self.restart_process = restart_process
        self.wtdata = wtdata
        self.uuid = wtdata["uuid"]
        self.name = wtdata["name"]

    def save(self):
        """Save Nodetree.

        - Update uuid for links. Build compressed nodes for nodetree.
        - Analysis connectivity
        - Check exist in database or not. If not in database, save directly.
        - If in database, analyze the difference, save accordingly.
        """
        self.build_node_link()
        self.build_connectivity()
        if self.exist_in_db() or self.restart_process is not None:
            new_nodes, modified_nodes, update_metadata = self.check_diff(
                self.restart_process
            )
            print("modified_nodes: {}".format(modified_nodes))
            self.reset_nodes(modified_nodes)
        self.insert_nodetree_to_db()

    def build_node_link(self):
        """Create links for nodes.
        Create the links for node inputs using:
        1) nodetree links

        """
        # reset node input links
        for name, node in self.wtdata["nodes"].items():
            for input in node["inputs"]:
                input["links"] = []
            for output in node["outputs"]:
                output["links"] = []
        for link in self.wtdata["links"]:
            to_socket = [
                socket
                for socket in self.wtdata["nodes"][link["to_node"]]["inputs"]
                if socket["name"] == link["to_socket"]
            ][0]
            from_socket = [
                socket
                for socket in self.wtdata["nodes"][link["from_node"]]["outputs"]
                if socket["name"] == link["from_socket"]
            ][0]
            to_socket["links"].append(link)
            from_socket["links"].append(link)

    def insert_nodetree_to_db(self):
        """Save a new nodetree in the database.

        - nodetree
        - all nodes
        """
        from aiida.orm.utils.serialize import serialize

        # pprint(self.wtdata)
        self.wtdata["created"] = datetime.datetime.utcnow()
        self.wtdata["lastUpdate"] = datetime.datetime.utcnow()
        self.process.base.extras.set("worktree", serialize(self.wtdata))

    def reset_nodes(self, nodes):
        """Reset nodes

        Args:
            nodes (list): a list of node names.
        """

        for name in nodes:
            self.append_message_to_queue(
                f"node,{name}:RESET",
            )

    def append_message_to_queue(self, message):
        queue = self.process.base.extras.get("worktree_queue", [])
        queue.append(message)
        self.process.base.extras.set("worktree_queue", queue)

    def set_nodes_action(self, action):
        """Set node action."""
        for name, node in self.wtdata["nodes"].items():
            # print("Reset node: {}".format(node))
            node["action"] = action

    def get_wtdata_from_db(self, process=None):
        from aiida.orm.utils.serialize import deserialize_unsafe

        process = self.process if process is None else process
        wtdata = process.base.extras.get("worktree", None)
        if wtdata is None:
            print("No worktree data found in the process node.")
            return
        wtdata = deserialize_unsafe(wtdata)
        return wtdata

    def check_diff(self, restart_process=None):
        """Find difference between nodetree and its database.

        Returns:
            new_nodes: new nodes
            modified_nodes: modified nodes
        """
        from node_graph.analysis import DifferenceAnalysis

        wt1 = self.get_wtdata_from_db(restart_process)
        dc = DifferenceAnalysis(nt1=wt1, nt2=self.wtdata)
        (
            new_nodes,
            modified_nodes,
            update_metadata,
        ) = dc.build_difference()
        return new_nodes, modified_nodes, update_metadata

    def exist_in_db(self):
        """Check nodetree exist in database or not.

        Returns:
            bool: _description_
        """
        if self.process.base.extras.get("worktree", None) is not None:
            return True
        return False

    def build_connectivity(self):
        """Analyze the connectivity of nodetree and save it into dict."""
        from node_graph.analysis import ConnectivityAnalysis

        nc = ConnectivityAnalysis(self.wtdata)
        self.wtdata["connectivity"] = nc.build_connectivity()
