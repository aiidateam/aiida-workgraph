from typing import Optional, Dict, Tuple, List

# import datetime
from aiida.orm import ProcessNode
from aiida.orm.utils.serialize import serialize, deserialize_unsafe


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
        for name, task in self.wgdata["tasks"].items():
            for wait_task in task["wait"]:
                if wait_task in self.wgdata["tasks"]:
                    self.wgdata["links"].append(
                        {
                            "from_node": wait_task,
                            "from_socket": "_wait",
                            "to_node": name,
                            "to_socket": "_wait",
                        }
                    )

    def clean_hanging_links(self) -> None:
        """Clean hanging links in the workgraph."""
        for link in self.wgdata["links"][:]:  # Iterate over a shallow copy of the list
            if (
                link["from_node"] not in self.wgdata["tasks"]
                or link["to_node"] not in self.wgdata["tasks"]
            ):
                self.wgdata["links"].remove(link)

    def save(self) -> None:
        """Save workgraph.

        - Update uuid for links. Build compressed tasks for workgraph.
        - Analysis connectivity
        - Check exist in database or not. If not in database, save directly.
        - If in database, analyze the difference, save accordingly.
        """
        self.build_task_link()
        self.build_connectivity()
        if self.exist_in_db() or self.restart_process is not None:
            new_tasks, modified_tasks, update_metadata = self.check_diff(
                self.restart_process
            )
            print("modified_tasks: {}".format(modified_tasks))
            self.reset_tasks(modified_tasks)
        self.insert_workgraph_to_db()

    def build_task_link(self) -> None:
        """Create links for tasks.
        Create the links for task inputs using:
        1) workgraph links

        """
        # reset task input links
        for name, task in self.wgdata["tasks"].items():
            for input in task["inputs"]:
                input["links"] = []
            for output in task["outputs"]:
                output["links"] = []
        for link in self.wgdata["links"]:
            to_socket = [
                socket
                for socket in self.wgdata["tasks"][link["to_node"]]["inputs"]
                if socket["name"] == link["to_socket"]
            ][0]
            from_socket = [
                socket
                for socket in self.wgdata["tasks"][link["from_node"]]["outputs"]
                if socket["name"] == link["from_socket"]
            ][0]
            to_socket["links"].append(link)
            from_socket["links"].append(link)

    def insert_workgraph_to_db(self) -> None:
        """Save a new workgraph in the database.

        - workgraph
        - all tasks
        """
        from aiida_workgraph.utils import workgraph_to_short_json

        # pprint(self.wgdata)
        # self.wgdata["created"] = datetime.datetime.utcnow()
        # self.wgdata["lastUpdate"] = datetime.datetime.utcnow()
        short_wgdata = workgraph_to_short_json(self.wgdata)
        self.process.base.extras.set("_workgraph_short", short_wgdata)
        self.save_task_states()
        for name, task in self.wgdata["tasks"].items():
            self.wgdata["tasks"][name] = serialize(task)
        # nodes is a copy of tasks, so we need to pop it out
        self.wgdata.pop("nodes")
        self.wgdata["error_handlers"] = serialize(self.wgdata["error_handlers"])
        self.process.base.extras.set("_workgraph", self.wgdata)

    def save_task_states(self) -> Dict:
        """Get task states."""

        task_states = {}
        task_processes = {}
        task_actions = {}
        for name, task in self.wgdata["tasks"].items():
            task_states[f"_task_state_{name}"] = task["state"]
            task_processes[f"_task_process_{name}"] = task["process"]
            task_actions[f"_task_action_{name}"] = task["action"]
        self.process.base.extras.set_many(task_states)
        self.process.base.extras.set_many(task_processes)
        self.process.base.extras.set_many(task_actions)

        return task_states

    def reset_tasks(self, tasks: List[str]) -> None:
        """Reset tasks

        Args:
            tasks (list): a list of task names.
        """
        from aiida_workgraph.utils.control import create_task_action

        # print("process state: ", self.process.process_state.value.upper())
        if (
            self.process.process_state is None
            or self.process.process_state.value.upper() == "CREATED"
        ):
            for name in tasks:
                self.wgdata["tasks"][name]["state"] = "PLANNED"
                self.wgdata["tasks"][name]["process"] = serialize(None)
                self.wgdata["tasks"][name]["result"] = None
                names = self.wgdata["connectivity"]["child_node"][name]
                for name in names:
                    self.wgdata["tasks"][name]["state"] = "PLANNED"
                    self.wgdata["tasks"][name]["result"] = None
                    self.wgdata["tasks"][name]["process"] = serialize(None)
        else:
            create_task_action(self.process.pk, tasks=tasks, action="reset")

    def set_tasks_action(self, action: str) -> None:
        """Set task action."""
        for name, task in self.wgdata["tasks"].items():
            # print("Reset task: {}".format(task))
            task["action"] = action

    def get_wgdata_from_db(
        self, process: Optional[ProcessNode] = None
    ) -> Optional[Dict]:

        process = self.process if process is None else process
        wgdata = process.base.extras.get("_workgraph", None)
        if wgdata is None:
            print("No workgraph data found in the process node.")
            return
        for name, task in wgdata["tasks"].items():
            wgdata["tasks"][name] = deserialize_unsafe(task)
        # also make a alias for nodes
        wgdata["nodes"] = wgdata["tasks"]
        wgdata["error_handlers"] = deserialize_unsafe(wgdata["error_handlers"])
        return wgdata

    def check_diff(
        self, restart_process: Optional[ProcessNode] = None
    ) -> Tuple[List[str], List[str], Dict]:
        """Find difference between workgraph and its database.

        Returns:
            new_tasks: new tasks
            modified_tasks: modified tasks
        """
        from node_graph.analysis import DifferenceAnalysis

        wg1 = self.get_wgdata_from_db(restart_process)
        dc = DifferenceAnalysis(nt1=wg1, nt2=self.wgdata)
        (
            new_tasks,
            modified_tasks,
            update_metadata,
        ) = dc.build_difference()
        return new_tasks, modified_tasks, update_metadata

    def exist_in_db(self) -> bool:
        """Check workgraph exist in database or not.

        Returns:
            bool: _description_
        """
        if self.process.base.extras.get("_workgraph", None) is not None:
            return True
        return False

    def build_connectivity(self) -> None:
        """Analyze the connectivity of workgraph and save it into dict."""
        from node_graph.analysis import ConnectivityAnalysis

        self.wgdata["nodes"] = self.wgdata["tasks"]
        nc = ConnectivityAnalysis(self.wgdata)
        self.wgdata["connectivity"] = nc.build_connectivity()
