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

    def save(self, update_state: bool = True) -> None:
        """Save workgraph."""

        self.separate_workgraph_data()
        if update_state:
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
