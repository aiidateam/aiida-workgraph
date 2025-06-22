from typing import Optional, Dict
from aiida.orm import ProcessNode, load_node
from aiida.orm.utils.serialize import serialize


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
        self.restart_process = (
            load_node(restart_process)
            if isinstance(restart_process, int)
            else restart_process
        )
        wgdata.setdefault("uuid", "")
        wgdata.setdefault("tasks", {})
        wgdata.setdefault("links", [])
        wgdata.setdefault("error_handlers", {})
        wgdata.setdefault("meta_sockets", {})
        self.wgdata = wgdata
        self.name = wgdata["name"]
        self.clean_hanging_links()

    def clean_hanging_links(self) -> None:
        """Clean hanging links in the workgraph."""
        for link in self.wgdata["links"][:]:  # Iterate over a shallow copy of the list
            if (
                link["from_node"] not in self.wgdata["tasks"]
                or link["to_node"] not in self.wgdata["tasks"]
            ):
                self.wgdata["links"].remove(link)

    def analyze(self) -> None:
        """Analyze workgraph.

        - Update uuid for links. Build compressed tasks for workgraph.
        - Analysis connectivity
        - Check exist in database or not. If not in database, save directly.
        - If in database, analyze the difference, save accordingly.
        """
        self.build_task_link()
        self.assign_zone()
        self.update_parent_task()
        self.find_all_zones_inputs()

    def save(self) -> None:
        """Save workgraph."""
        self.analyze()
        self.insert_workgraph_to_db()

    def build_task_link(self) -> None:
        """Create links for tasks.
        Create the links for task inputs using:
        1) workgraph links

        """
        # create a `input_links` to store the input links for each task
        for task in self.wgdata["tasks"].values():
            task["input_links"] = {}
        for link in self.wgdata["links"]:
            task = self.wgdata["tasks"][link["to_node"]]
            if link["to_socket"] not in task["input_links"]:
                task["input_links"][link["to_socket"]] = []
            task["input_links"][link["to_socket"]].append(link)

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
        for _, links in self.wgdata["tasks"][name]["input_links"].items():
            for link in links:
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
                for _, links in self.wgdata["tasks"][child_task]["input_links"].items():
                    for link in links:
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
        self.separate_workgraph_data()
        self.process.task_states = self.task_states
        self.process.task_processes = self.task_processes
        self.process.task_actions = self.task_actions
        self.process.task_executors = self.task_executors
        self.process.task_error_handlers = self.task_error_handlers
        self.process.workgraph_data = self.wgdata
        self.process.workgraph_data_short = self.short_wgdata
        self.process.workgraph_error_handlers = self.workgraph_error_handlers

    def separate_workgraph_data(self) -> None:
        """Separate the data into different parts and store them in the process node.
        WorkGraph data:
        - workgraph_data
        - task_executors, which could contain pickled binary data
        - task_error_handlers, which could contain pickled binary data
        - workgraph_error_handlers, which could contain pickled binary data

        Runtime information:
        - task_states
        - task_processes
        - task_actions

        """
        from aiida_workgraph.utils import (
            workgraph_to_short_json,
            serialize_input_values_recursively,
        )

        self.task_states = {}
        self.task_processes = {}
        self.task_actions = {}
        self.task_executors = {}
        self.task_error_handlers = {}
        self.workgraph_error_handlers = {}
        self.short_wgdata = workgraph_to_short_json(self.wgdata)
        for name, task in self.wgdata["tasks"].items():
            self.task_states[name] = task["state"]
            self.task_processes[name] = task["process"]
            self.task_actions[name] = task["action"]
            self.task_executors[name] = task.pop("executor", None)
            self.task_error_handlers[name] = task.pop("error_handlers", {})
            serialize_input_values_recursively(task["inputs"])
        self.workgraph_error_handlers = self.wgdata.pop("error_handlers")
        self.wgdata["meta_sockets"] = serialize(self.wgdata["meta_sockets"])
