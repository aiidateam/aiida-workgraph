from typing import Optional, Dict, Tuple, List
import datetime
from aiida.orm import ProcessNode


class WorkGraphSaver:
    """Save a workgraph to the database."""

    def __init__(
        self,
        process: ProcessNode,
        wgdata: Dict,
        restart_process: Optional[ProcessNode] = None,
    ) -> None:
        """Init WorkGraphSaver.

        Args:
            wgdata (dict): data of workgraph to be launched.
        """
        self.process = process
        self.restart_process = restart_process
        self.wgdata = wgdata
        self.uuid = wgdata["uuid"]
        self.name = wgdata["name"]
        self.wait_to_link()
        self.clean_hanging_links()

    def wait_to_link(self) -> None:
        """Convert wait attribute to link."""
        for name, node in self.wgdata["nodes"].items():
            for wait_node in node["wait"]:
                if wait_node in self.wgdata["nodes"]:
                    self.wgdata["links"].append(
                        {
                            "from_node": wait_node,
                            "from_socket": "_wait",
                            "to_node": name,
                            "to_socket": "_wait",
                        }
                    )

    def clean_hanging_links(self) -> None:
        """Clean hanging links in the workgraph."""
        for link in self.wgdata["links"][:]:  # Iterate over a shallow copy of the list
            if (
                link["from_node"] not in self.wgdata["nodes"]
                or link["to_node"] not in self.wgdata["nodes"]
            ):
                self.wgdata["links"].remove(link)

    def save(self) -> None:
        """Save workgraph.

        - Update uuid for links. Build compressed nodes for workgraph.
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
        self.insert_workgraph_to_db()

    def build_node_link(self) -> None:
        """Create links for nodes.
        Create the links for node inputs using:
        1) workgraph links

        """
        # reset node input links
        for name, node in self.wgdata["nodes"].items():
            for input in node["inputs"]:
                input["links"] = []
            for output in node["outputs"]:
                output["links"] = []
        for link in self.wgdata["links"]:
            to_socket = [
                socket
                for socket in self.wgdata["nodes"][link["to_node"]]["inputs"]
                if socket["name"] == link["to_socket"]
            ][0]
            from_socket = [
                socket
                for socket in self.wgdata["nodes"][link["from_node"]]["outputs"]
                if socket["name"] == link["from_socket"]
            ][0]
            to_socket["links"].append(link)
            from_socket["links"].append(link)

    def insert_workgraph_to_db(self) -> None:
        """Save a new workgraph in the database.

        - workgraph
        - all nodes
        """
        from aiida.orm.utils.serialize import serialize

        # pprint(self.wgdata)
        self.wgdata["created"] = datetime.datetime.utcnow()
        self.wgdata["lastUpdate"] = datetime.datetime.utcnow()
        self.process.base.extras.set("workgraph", serialize(self.wgdata))

    def reset_nodes(self, nodes: List[str]) -> None:
        """Reset nodes

        Args:
            nodes (list): a list of node names.
        """

        for name in nodes:
            self.append_message_to_queue(
                f"node,{name}:RESET",
            )

    def append_message_to_queue(self, message: str) -> None:
        queue = self.process.base.extras.get("workgraph_queue", [])
        queue.append(message)
        self.process.base.extras.set("workgraph_queue", queue)

    def set_nodes_action(self, action: str) -> None:
        """Set node action."""
        for name, node in self.wgdata["nodes"].items():
            # print("Reset node: {}".format(node))
            node["action"] = action

    def get_wgdata_from_db(
        self, process: Optional[ProcessNode] = None
    ) -> Optional[Dict]:
        from aiida.orm.utils.serialize import deserialize_unsafe

        process = self.process if process is None else process
        wgdata = process.base.extras.get("workgraph", None)
        if wgdata is None:
            print("No workgraph data found in the process node.")
            return
        wgdata = deserialize_unsafe(wgdata)
        return wgdata

    def check_diff(
        self, restart_process: Optional[ProcessNode] = None
    ) -> Tuple[List[str], List[str], Dict]:
        """Find difference between workgraph and its database.

        Returns:
            new_nodes: new nodes
            modified_nodes: modified nodes
        """
        from node_graph.analysis import DifferenceAnalysis

        wg1 = self.get_wgdata_from_db(restart_process)
        dc = DifferenceAnalysis(nt1=wg1, nt2=self.wgdata)
        (
            new_nodes,
            modified_nodes,
            update_metadata,
        ) = dc.build_difference()
        return new_nodes, modified_nodes, update_metadata

    def exist_in_db(self) -> bool:
        """Check workgraph exist in database or not.

        Returns:
            bool: _description_
        """
        if self.process.base.extras.get("workgraph", None) is not None:
            return True
        return False

    def build_connectivity(self) -> None:
        """Analyze the connectivity of workgraph and save it into dict."""
        from node_graph.analysis import ConnectivityAnalysis

        nc = ConnectivityAnalysis(self.wgdata)
        self.wgdata["connectivity"] = nc.build_connectivity()
