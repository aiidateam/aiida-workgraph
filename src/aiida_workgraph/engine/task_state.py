from __future__ import annotations
from typing import Optional, Tuple, List, Any
from aiida.orm.utils.serialize import serialize
from aiida_workgraph.orm.utils import deserialize_safe
from aiida.orm import ProcessNode, Data


class TaskStateManager:
    """
    Handles all low-level operations on tasks' states, runtime info,
    and relationships (parent/child).
    """

    def __init__(self, ctx_manager, logger, process, awaitable_manager):
        """
        :param ctx_manager: Context manager holding `ctx` (containing tasks, connectivity, etc).
        :param logger: Logger instance.
        :param process: The current AiiDA process.
        :param awaitable_manager: Manager that orchestrates async tasks/futures.
        """
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger
        self.process = process
        self.awaitable_manager = awaitable_manager

    # -------------------------
    # Basic get/set operations
    # -------------------------
    def get_task_runtime_info(self, name: str, key: str) -> Any:
        """Fetch a task runtime property (e.g. process, state, action)."""
        if key == "process":
            value = self.process.node.get_task_process(name)
            return deserialize_safe(value) if value else None
        elif key == "state":
            return self.process.node.get_task_state(name)
        elif key == "action":
            return self.process.node.get_task_action(name)
        elif key == "execution_count":
            return self.process.node.get_task_execution_count(name)
        else:
            raise ValueError(f"Invalid key: {key}")
        return None

    def set_task_runtime_info(self, name: str, key: str, value: Any) -> None:
        """Set a task runtime property (e.g. process, state, action).
        All the runtime info are store into the process node, which allow us
        access this info outside the engine
        """
        if key == "process":
            serialized = serialize(value)
            self.process.node.set_task_process(name, serialized)
        elif key == "state":
            self.process.node.set_task_state(name, value)
        elif key == "action":
            self.process.node.set_task_action(name, value)
        elif key == "execution_count":
            self.process.node.set_task_execution_count(name, value)
        elif key == "map_info":
            self.process.node.set_task_map_info(name, value)
        else:
            raise ValueError(f"Invalid key: {key}")

    def set_tasks_state(self, tasks: List[str], value: str) -> None:
        """
        Set the state for a list of tasks (and their children) to `value`.
        Typically used for skip or reset tasks.
        """
        for name in tasks:
            self.set_task_runtime_info(name, "state", value)
            if hasattr(self.process.wg.tasks[name], "children"):
                self.set_tasks_state(
                    [task.name for task in self.process.wg.tasks[name].children], value
                )
            # TODO should we also reset the mapped tasks?

    def update_task_state(self, name: str, success=True) -> None:
        """Update task state when the task is finished."""
        task = self.process.wg.tasks[name]
        self.ctx._task_results.setdefault(name, {})
        if success:
            node = self.get_task_runtime_info(name, "process")
            if isinstance(node, ProcessNode):
                state = node.process_state.value.upper()
                if node.is_finished_ok:
                    self.set_task_runtime_info(task.name, "state", state)
                    if task.node_type.upper() == "WORKGRAPH":
                        # expose the outputs of all the tasks in the workgraph
                        outgoing = node.base.links.get_outgoing()
                        for link in outgoing.all():
                            if isinstance(link.node, ProcessNode) and getattr(
                                link.node, "process_state", False
                            ):
                                self.ctx._task_results[name][
                                    link.link_label
                                ] = link.node.outputs
                    else:
                        self.ctx._task_results[name] = node.outputs
                        # self.ctx._new_data[name] = self.ctx._task_results[name]
                    self.set_task_runtime_info(task.name, "state", "FINISHED")
                    self.update_meta_tasks(name)
                    self.process.report(
                        f"Task: {name}, type: {task.node_type}, finished."
                    )
                # all other states are considered as failed
                else:
                    self.ctx._task_results[name] = node.outputs
                    self.on_task_failed(name)
            elif isinstance(node, Data):
                #
                output_name = [
                    output_name
                    for output_name in task.outputs._get_keys()
                    if output_name not in ["_wait", "_outputs"]
                ][0]
                self.ctx._task_results[name] = {output_name: node}
                self.set_task_runtime_info(task.name, "state", "FINISHED")
                self.update_meta_tasks(name)
                self.process.report(f"Task: {name} finished.")
        else:
            self.on_task_failed(name)
        # After finishing, inform the parent
        self.update_parent_task_state(name)

    def update_normal_task_state(self, name, results, success=True):
        """Set the results of a normal task.
        A normal task is created by decorating a function with @task().
        """
        from aiida_workgraph.config import builtin_outputs

        builtin_output_names = [output["name"] for output in builtin_outputs]

        if success:
            task = self.process.wg.tasks[name]
            if isinstance(results, tuple):
                # there are two built-in outputs: _wait and _outputs
                if len(task.outputs) - 2 != len(results):
                    self.on_task_failed(name)
                    return self.process.exit_codes.OUTPUS_NOT_MATCH_RESULTS
                output_names = [
                    output._name
                    for output in task.outputs
                    if output._name not in builtin_output_names
                ]
                for i, output_name in enumerate(output_names):
                    self.ctx._task_results[name][output_name] = results[i]
            elif isinstance(results, dict):
                self.ctx._task_results[name] = results
            else:
                output_names = [
                    output_name
                    for output_name in task.outputs._get_keys()
                    if output_name not in ["_wait", "_outputs"]
                ]
                # some task does not have any output
                if len(output_names) == 1:
                    self.ctx._task_results[name][output_names[0]] = results
                    if isinstance(results, Data):
                        results.store()
                        self.set_task_runtime_info(task.name, "process", results)
                elif len(output_names) > 1:
                    self.process.exit_codes.OUTPUS_NOT_MATCH_RESULTS
            self.update_meta_tasks(name)
            self.set_task_runtime_info(name, "state", "FINISHED")
            self.process.report(f"Task: {name} finished.")
        else:
            self.on_task_failed(name)
        self.update_parent_task_state(name)

    def update_meta_tasks(self, name: str) -> None:
        """Export task results to the context based on context mapping."""
        from aiida_workgraph.utils import update_nested_dict, get_nested_dict

        for link in self.process.wg.links:
            if link.from_node.name == name and link.to_node.name in [
                "graph_ctx",
                "graph_outputs",
            ]:
                key = link.to_socket._scoped_name
                result_key = link.from_socket._scoped_name
                result = get_nested_dict(
                    self.ctx._task_results[name], result_key, default=None
                )
                update_nested_dict(
                    self.ctx._task_results[link.to_node.name], key, result
                )

    # --------------------------------------------------
    # Reset & removing from executed tasks
    # --------------------------------------------------
    def reset_task(
        self,
        name: str,
        reset_process: bool = True,
        recursive: bool = True,
        reset_execution_count: bool = True,
    ) -> None:
        """
        Reset the task's state to PLANNED, optionally clearing the process reference
        and recursing to children. If the task is a WHILE, reset its execution_count.
        """
        self.logger.debug(f"Resetting task {name}.")
        self.set_task_runtime_info(name, "state", "PLANNED")
        if reset_process:
            self.set_task_runtime_info(name, "process", None)
        self.remove_executed_task(name)

        node_type = self.process.wg.tasks[name].node_type.upper()
        if node_type == "WHILE":
            if reset_execution_count:
                self.set_task_runtime_info(name, "execution_count", 0)
            for child_task in self.process.wg.tasks[name].children:
                self.reset_task(child_task.name, reset_process=False, recursive=False)
        elif node_type in ["IF", "ZONE"]:
            for child_task in self.process.wg.tasks[name].children:
                self.reset_task(child_task.name, reset_process=False, recursive=False)

        if recursive:
            # reset its child tasks
            child_names = self.process.wg.connectivity["child_node"][name]
            for child_name in child_names:
                self.reset_task(child_name, recursive=False)

        self.logger.debug(f"Task {name} was reset.")

    def remove_executed_task(self, name: str) -> None:
        """
        Remove tasks from `ctx._executed_tasks` if they match this name (or name.*).
        """
        self.ctx._executed_tasks = [
            label for label in self.ctx._executed_tasks if label.split(".")[0] != name
        ]

    # --------------------------------------------------
    # Checking readiness, finishing, failures, etc.
    # --------------------------------------------------
    def is_task_ready_to_run(self, name: str) -> Tuple[bool, Optional[str]]:
        """
        Check if the task is ready to run. We consider parent states, input tasks, etc.
        For tasks inside a ZONE or with a parent task, we require the parent
        to be in a running state, and the zone's input tasks finished or failed.
        """
        parent_task = self.process.wg.tasks[name].parent_task
        parent_states = [True, True]

        # If the task has a parent zone
        if parent_task:
            state = self.get_task_runtime_info(parent_task.name, "state")
            if state not in ["RUNNING"]:
                parent_states[1] = False

        # Check input tasks from the zone connectivity
        for child_task_name in self.process.wg.connectivity["zone"][name][
            "input_tasks"
        ]:
            child_state = self.get_task_runtime_info(child_task_name, "state")
            if child_state not in ["FINISHED", "SKIPPED", "FAILED"]:
                parent_states[0] = False
                break
        return all(parent_states), parent_states

    # --------------------------------------------------
    # Relationship & parent updates
    # --------------------------------------------------
    def on_task_failed(self, name: str) -> None:
        """
        Mark a task as FAILED, skip its children, and run any error handlers.
        """
        task_type = self.process.wg.tasks[name].node_type
        self.set_task_runtime_info(name, "state", "FAILED")
        self.set_tasks_state(
            self.process.wg.connectivity["child_node"][name], "SKIPPED"
        )
        msg = f"Task, {name}, type: {task_type}, failed."
        process = self.get_task_runtime_info(name, "process")
        if isinstance(process, ProcessNode):
            msg += f" Error message: {process.exit_message}"
        self.process.report(msg)
        self.process.error_handler_manager.run_error_handlers(name)

    def update_parent_task_state(self, name: str) -> None:
        """
        If a task has a parent (WHILE, IF, ZONE, MAP), notify the parent to update
        its own state. Also handle mapped tasks referencing a 'map_data.parent' node.
        """
        parent_task = self.process.wg.tasks[name].parent_task
        if parent_task:
            node_type = parent_task.node_type.upper()
            if node_type == "WHILE":
                self.update_while_task_state(parent_task.name)
            elif node_type in ["IF", "ZONE", "MAP"]:
                self.update_zone_task_state(parent_task.name)

        # If the task is a mapped child, update its parent's "template" (the original map node)
        if self.process.wg.tasks[name].map_data:
            map_parent = self.process.wg.tasks[name].map_data["parent"]
            self.update_template_task_state(map_parent)

    def update_while_task_state(self, name: str) -> None:
        """
        Called when a child of a WHILE task finishes. If all children are done, we decide
        whether to reset for the next iteration or finalize the WHILE.
        """
        finished, _ = self.are_childen_finished(name)

        if finished:
            self.process.report(
                f"Wihle Task {name}: this iteration finished. Try to reset for the next iteration."
            )
            # reset the condition tasks
            for link in self.process.wg.tasks[name].inputs.conditions._links:
                self.reset_task(link.from_node.name, recursive=False)
            # reset the task and all its children, so that the task can run again
            # do not reset the execution count
            self.reset_task(name, reset_execution_count=False)

    def update_zone_task_state(self, name: str) -> None:
        """
        Update the state of an IF or ZONE block. Mark it FINISHED if children are done.
        """
        finished, _ = self.are_childen_finished(name)
        if finished:
            self.set_task_runtime_info(name, "state", "FINISHED")
            self.process.report(f"Task: {name} finished.")
            self.update_parent_task_state(name)

    def update_template_task_state(self, name: str) -> None:
        """Update the template task state.
        1) check if all child tasks are finished.
        2) gather the results of all the mapped tasks.
        3) update the parent task state.
        """
        finished, _ = self.are_childen_finished(name)
        if finished:
            # gather the results of all the mapped tasks
            results = {}
            for prefix, mapped_task in self.process.wg.tasks[name].mapped_tasks.items():
                for output in mapped_task.outputs:
                    if output._name in self.ctx._task_results[mapped_task.name]:
                        results.setdefault(output._name, {})
                        results[output._name][prefix] = self.ctx._task_results[
                            mapped_task.name
                        ][output._name]
            self.ctx._task_results[name] = results
            self.set_task_runtime_info(name, "state", "FINISHED")
            self.update_meta_tasks(name)
            self.process.report(f"Task: {name} finished.")
            self.update_parent_task_state(name)

    def are_childen_finished(self, name: str) -> tuple[bool, Any]:
        """Check if the child tasks are finished."""
        task = self.process.wg.tasks[name]
        finished = True
        if hasattr(task, "children"):
            for child in task.children:
                if self.get_task_runtime_info(child.name, "state") not in [
                    "FINISHED",
                    "SKIPPED",
                    "FAILED",
                ]:
                    finished = False
                    break
        # check the mapped tasks
        mapped_tasks = task.mapped_tasks or {}
        for mapped_task in mapped_tasks.values():
            if self.get_task_runtime_info(mapped_task.name, "state") not in [
                "FINISHED",
                "SKIPPED",
                "FAILED",
            ]:
                finished = False
                break
        return finished, None
