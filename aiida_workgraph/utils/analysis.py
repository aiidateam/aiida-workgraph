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
                else:
                    raise ValueError(
                        "Task {} wait for a non-exist task {}".format(name, wait_task)
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
        self.assign_zone()
        self.build_connectivity()
        self.update_parent_task()
        self.find_all_zones_inputs()
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
            for _, input in task["inputs"].items():
                input["links"] = []
            for _, output in task["outputs"].items():
                output["links"] = []
        for link in self.wgdata["links"]:
            to_socket = [
                socket
                for name, socket in self.wgdata["tasks"][link["to_node"]][
                    "inputs"
                ].items()
                if name == link["to_socket"]
            ][0]
            from_socket = [
                socket
                for name, socket in self.wgdata["tasks"][link["from_node"]][
                    "outputs"
                ].items()
                if name == link["from_socket"]
            ][0]
            to_socket["links"].append(link)
            from_socket["links"].append(link)

    def assign_zone(self) -> None:
        """Assign zone for each task."""
        # assign parent_task for each task
        for name, task in self.wgdata["tasks"].items():
            for child_task in task["children"]:
                self.wgdata["tasks"][child_task]["parent_task"][0] = name

    def update_parent_task(self) -> None:
        """Recursively update the list of parent tasks for each task in wgdata."""

        def get_all_parents(task_name):
            """Recursively collect all parent tasks for a given task."""
            task = self.wgdata["tasks"][task_name]
            parent_names = []
            while task["parent_task"][0] is not None:
                parent_task_name = task["parent_task"][0]
                parent_names.append(parent_task_name)
                task = self.wgdata["tasks"][parent_task_name]
            parent_names.append(None)
            return parent_names

        # Loop through all tasks and update their parent task lists recursively
        for name, task in self.wgdata["tasks"].items():
            task["parent_task"] = get_all_parents(name)

    def find_all_zones_inputs(self) -> None:
        for name in self.wgdata["tasks"]:
            self.find_zone_inputs(name)

    def find_zone_inputs(self, name: str) -> None:
        """Find the input and outputs tasks for the zone."""
        task = self.wgdata["tasks"][name]
        input_tasks = []
        for _, input in self.wgdata["tasks"][name]["inputs"].items():
            for link in input["links"]:
                input_tasks.append(link["from_node"])
        # find all the input tasks
        for child_task in task["children"]:
            # if the child task is a zone
            if self.wgdata["tasks"][child_task]["children"]:
                # find the input tasks of the child task zone
                self.find_zone_inputs(child_task)
                # find all the input tasks which outside the while zone
                for child_task1 in self.wgdata["connectivity"]["zone"][child_task][
                    "input_tasks"
                ]:
                    if child_task1 not in task["children"]:
                        input_tasks.append(child_task1)
            else:
                # if the child task is not a zone, get the input tasks of the child task
                # find all the input tasks which outside the while zone
                for _, input in self.wgdata["tasks"][child_task]["inputs"].items():
                    for link in input["links"]:
                        input_tasks.append(link["from_node"])
        # find the input tasks which are not in the zone
        new_input_tasks = []
        for input_task in input_tasks:
            task_to_check = [input_task]
            task_to_check.extend(self.wgdata["tasks"][input_task]["parent_task"])
            not_in_the_zone = True
            for task_name in task_to_check:
                if task_name in task["children"]:
                    not_in_the_zone = False
                    break
            if not_in_the_zone:
                new_input_tasks.append(input_task)
        # find the parent task of the input tasks
        final_input_tasks = []
        parent_tasks = task["parent_task"]
        for input_task in new_input_tasks:
            # find the first parent task of this two task which are the same
            for parent_task in parent_tasks:
                if parent_task in self.wgdata["tasks"][input_task]["parent_task"]:
                    break
            # add the input task to the parent task
            if parent_task in self.wgdata["tasks"][input_task]["parent_task"]:
                index = self.wgdata["tasks"][input_task]["parent_task"].index(
                    parent_task
                )
            else:
                index = 0
            if index == 0:
                final_input_tasks.append(input_task)
            else:
                final_input_tasks.append(
                    self.wgdata["tasks"][input_task]["parent_task"][index - 1]
                )
        # remove the duplicate tasks
        final_input_tasks = list(set(final_input_tasks))
        self.wgdata["connectivity"]["zone"][name] = {"input_tasks": final_input_tasks}

    def insert_workgraph_to_db(self) -> None:
        """Save a new workgraph in the database.

        - workgraph
        - all tasks
        """
        from aiida_workgraph.utils import workgraph_to_short_json
        import inspect
        from aiida_workgraph.orm.function_data import PickledLocalFunction

        # pprint(self.wgdata)
        # self.wgdata["created"] = datetime.datetime.utcnow()
        # self.wgdata["lastUpdate"] = datetime.datetime.utcnow()
        short_wgdata = workgraph_to_short_json(self.wgdata)
        self.process.base.extras.set("_workgraph_short", short_wgdata)
        self.save_task_states()
        for name, task in self.wgdata["tasks"].items():
            for _, input in task["inputs"].items():
                if input["property"] is None:
                    continue
                prop = input["property"]
                if inspect.isfunction(prop["value"]):
                    prop["value"] = PickledLocalFunction(prop["value"]).store()
            self.wgdata["tasks"][name] = serialize(task)
        # nodes is a copy of tasks, so we need to pop it out
        self.wgdata["error_handlers"] = serialize(self.wgdata["error_handlers"])
        self.wgdata["context"] = serialize(self.wgdata["context"])
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
        # change tasks to nodes for DifferenceAnalysis
        wg1["nodes"] = wg1.pop("tasks")
        self.wgdata["nodes"] = self.wgdata.pop("tasks")
        dc = DifferenceAnalysis(nt1=wg1, nt2=self.wgdata)
        (
            new_tasks,
            modified_tasks,
            update_metadata,
        ) = dc.build_difference()
        # change nodes back to tasks
        wg1["tasks"] = wg1.pop("nodes")
        self.wgdata["tasks"] = self.wgdata.pop("nodes")
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

        # ConnectivityAnalysis use nodes instead of tasks
        self.wgdata["nodes"] = self.wgdata.pop("tasks")
        nc = ConnectivityAnalysis(self.wgdata)
        self.wgdata["connectivity"] = nc.build_connectivity()
        self.wgdata["connectivity"]["zone"] = {}
        # change nodes back to tasks
        self.wgdata["tasks"] = self.wgdata.pop("nodes")
