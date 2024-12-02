from __future__ import annotations

from aiida.orm import ProcessNode, Data
from aiida.engine import run_get_node
from typing import Any, Dict, List, Optional, Tuple, Callable, Union, Sequence
from aiida_workgraph.utils import create_and_pause_process
from aiida_workgraph.task import Task
from aiida_workgraph.utils import get_nested_dict
from aiida.orm.utils.serialize import deserialize_unsafe, serialize
import asyncio
from aiida.engine.processes.exit_code import ExitCode
from aiida_workgraph.executors.monitors import monitor

MAX_NUMBER_AWAITABLES_MSG = "The maximum number of subprocesses has been reached: {}. Cannot launch the job: {}."

process_task_types = [
    "CALCJOB",
    "WORKCHAIN",
    "GRAPH_BUILDER",
    "WORKGRAPH",
    "PYTHONJOB",
    "SHELLJOB",
]


class TaskManager:
    """Manages task execution, state updates, and error handling."""

    def __init__(self, ctx_manager, logger, runner, process, awaitable_manager):
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger
        self.runner = runner
        self.process = process
        self.awaitable_manager = awaitable_manager

    def get_task(self, name: str):
        """Get task from the context."""
        task = Task.from_dict(self.ctx._tasks[name])
        # update task results
        for output in task.outputs:
            output.value = get_nested_dict(
                self.ctx._tasks[name]["results"],
                output.name,
                default=output.value,
            )
        return task

    def reset_task(
        self,
        name: str,
        reset_process: bool = True,
        recursive: bool = True,
        reset_execution_count: bool = True,
    ) -> None:
        """Reset task state and remove it from the executed task.
        If recursive is True, reset its child tasks."""

        self.set_task_state_info(name, "state", "PLANNED")
        if reset_process:
            self.set_task_state_info(name, "process", None)
        self.remove_executed_task(name)
        # self.logger.debug(f"Task {name} action: RESET.")
        # if the task is a while task, reset its child tasks
        if self.ctx._tasks[name]["metadata"]["node_type"].upper() == "WHILE":
            if reset_execution_count:
                self.ctx._tasks[name]["execution_count"] = 0
            for child_task in self.ctx._tasks[name]["children"]:
                self.reset_task(child_task, reset_process=False, recursive=False)
        elif self.ctx._tasks[name]["metadata"]["node_type"].upper() in [
            "IF",
            "ZONE",
        ]:
            for child_task in self.ctx._tasks[name]["children"]:
                self.reset_task(child_task, reset_process=False, recursive=False)
        if recursive:
            # reset its child tasks
            names = self.ctx._connectivity["child_node"][name]
            for name in names:
                self.reset_task(name, recursive=False)

    def remove_executed_task(self, name: str) -> None:
        """Remove labels with name from executed tasks."""
        self.ctx._executed_tasks = [
            label for label in self.ctx._executed_tasks if label.split(".")[0] != name
        ]

    def is_task_ready_to_run(self, name: str) -> Tuple[bool, Optional[str]]:
        """Check if the task ready to run.
        For normal task and a zone task, we need to check its input tasks in the connectivity["zone"].
        For task inside a zone, we need to check if the zone (parent task) is ready.
        """
        parent_task = self.ctx._tasks[name]["parent_task"]
        # input_tasks, parent_task, conditions
        parent_states = [True, True]
        # if the task belongs to a parent zone
        if parent_task[0]:
            state = self.get_task_state_info(parent_task[0], "state")
            if state not in ["RUNNING"]:
                parent_states[1] = False
        # check the input tasks of the zone
        # check if the zone input tasks are ready
        for child_task_name in self.ctx._connectivity["zone"][name]["input_tasks"]:
            if self.get_task_state_info(child_task_name, "state") not in [
                "FINISHED",
                "SKIPPED",
                "FAILED",
            ]:
                parent_states[0] = False
                break

        return all(parent_states), parent_states

    def set_task_results(self) -> None:
        for name, task in self.ctx._tasks.items():
            if self.get_task_state_info(name, "action").upper() == "RESET":
                self.reset_task(task["name"])
            self.update_task_state(name)

    def task_set_context(self, name: str) -> None:
        """Export task results to the context based on context mapping."""
        from aiida_workgraph.utils import update_nested_dict

        items = self.ctx._tasks[name]["context_mapping"]
        for key, result_name in items.items():
            result = self.ctx._tasks[name]["results"][result_name]
            update_nested_dict(self.ctx, key, result)

    def get_task_state_info(self, name: str, key: str) -> str:
        """Get task state info from ctx."""

        value = self.ctx._tasks[name].get(key, None)
        if key == "process" and value is not None:
            value = deserialize_unsafe(value)
        return value

    def set_task_state_info(self, name: str, key: str, value: any) -> None:
        """Set task state info to ctx and base.extras.
        We task state to the base.extras, so that we can access outside the engine"""

        if key == "process":
            value = serialize(value)
            self.process.node.base.extras.set(f"_task_{key}_{name}", value)
        else:
            self.process.node.base.extras.set(f"_task_{key}_{name}", value)
        self.ctx._tasks[name][key] = value

    def set_tasks_state(
        self, tasks: Union[List[str], Sequence[str]], value: str
    ) -> None:
        """Set tasks state"""
        for name in tasks:
            self.set_task_state_info(name, "state", value)
            if "children" in self.ctx._tasks[name]:
                self.set_tasks_state(self.ctx._tasks[name]["children"], value)

    def is_workgraph_finished(self) -> bool:
        """Check if the workgraph is finished.
        For `while` workgraph, we need check its conditions"""
        is_finished = True
        failed_tasks = []
        for name, task in self.ctx._tasks.items():
            # self.update_task_state(name)
            if self.get_task_state_info(task["name"], "state") in [
                "RUNNING",
                "CREATED",
                "PLANNED",
                "READY",
            ]:
                is_finished = False
            elif self.get_task_state_info(task["name"], "state") == "FAILED":
                failed_tasks.append(name)
        if is_finished:
            if self.ctx._workgraph["workgraph_type"].upper() == "WHILE":
                should_run = self.check_while_conditions()
                is_finished = not should_run
            if self.ctx._workgraph["workgraph_type"].upper() == "FOR":
                should_run = self.check_for_conditions()
                is_finished = not should_run
        if is_finished and len(failed_tasks) > 0:
            message = f"WorkGraph finished, but tasks: {failed_tasks} failed. Thus all their child tasks are skipped."
            self.process.report(message)
            result = ExitCode(302, message)
        else:
            result = None
        return is_finished, result

    def continue_workgraph(self) -> None:
        self.process.report("Continue workgraph.")
        task_to_run = []
        for name, task in self.ctx._tasks.items():
            # update task state
            if (
                self.get_task_state_info(task["name"], "state")
                in [
                    "CREATED",
                    "RUNNING",
                    "FINISHED",
                    "FAILED",
                    "SKIPPED",
                ]
                or name in self.ctx._executed_tasks
            ):
                continue
            ready, _ = self.is_task_ready_to_run(name)
            if ready:
                task_to_run.append(name)
        #
        self.process.report("tasks ready to run: {}".format(",".join(task_to_run)))
        self.run_tasks(task_to_run)

    def run_tasks(self, names: List[str], continue_workgraph: bool = True) -> None:
        """Run tasks.
        Task type includes: Node, Data, CalcFunction, WorkFunction, CalcJob, WorkChain, GraphBuilder,
        WorkGraph, PythonJob, ShellJob, While, If, Zone, GetContext, SetContext, Normal.

        """
        from aiida_workgraph.utils import (
            get_executor,
            update_nested_dict_with_special_keys,
        )

        for name in names:
            # skip if the max number of awaitables is reached
            task = self.ctx._tasks[name]
            if task["metadata"]["node_type"].upper() in process_task_types:
                if len(self.process._awaitables) >= self.ctx._max_number_awaitables:
                    print(
                        MAX_NUMBER_AWAITABLES_MSG.format(
                            self.ctx._max_number_awaitables, name
                        )
                    )
                    continue
            # skip if the task is already executed
            # or if the task is in a skippped state
            if name in self.ctx._executed_tasks or self.get_task_state_info(
                name, "state"
            ) in ["SKIPPED"]:
                continue
            self.ctx._executed_tasks.append(name)
            print("-" * 60)

            self.process.report(
                f"Run task: {name}, type: {task['metadata']['node_type']}"
            )
            executor, _ = get_executor(task["executor"])
            args, kwargs, var_args, var_kwargs, args_dict = self.get_inputs(name)
            for i, key in enumerate(self.ctx._tasks[name]["args"]):
                kwargs[key] = args[i]
            # update the port namespace
            kwargs = update_nested_dict_with_special_keys(kwargs)
            # kwargs["meta.label"] = name
            # output must be a Data type or a mapping of {string: Data}
            task["results"] = {}
            task_type = task["metadata"]["node_type"].upper()
            if task_type == "NODE":
                self.execute_node_task(
                    name, executor, kwargs, var_args, var_kwargs, continue_workgraph
                )
            elif task_type == "DATA":
                self.execute_data_task(name, executor, args, kwargs, continue_workgraph)
            elif task_type in ["CALCFUNCTION", "WORKFUNCTION"]:
                self.execute_function_task(
                    name, executor, kwargs, var_kwargs, continue_workgraph
                )
            elif task_type in ["CALCJOB", "WORKCHAIN"]:
                self.execute_process_task(name, executor, kwargs)
            elif task_type == "PYTHONJOB":
                self.execute_python_job_task(task, kwargs, var_kwargs)
            elif task_type == "GRAPH_BUILDER":
                self.execute_graph_builder_task(
                    task, executor, kwargs, var_args, var_kwargs
                )
            elif task_type in ["WORKGRAPH"]:
                self.execute_workgraph_task(task, kwargs)
            elif task_type == "SHELLJOB":
                self.execute_shell_job_task(task, kwargs)
            elif task_type == "WHILE":
                self.execute_while_task(task)
            elif task_type == "IF":
                self.execute_if_task(task)
            elif task_type == "ZONE":
                self.execute_zone_task(task)
            elif task_type == "GET_CONTEXT":
                self.execute_get_context_task(task, kwargs)
            elif task_type == "SET_CONTEXT":
                self.execute_set_context_task(task, kwargs)
            elif task_type == "AWAITABLE":
                self.execute_awaitable_task(
                    task, executor, args, kwargs, var_args, var_kwargs
                )
            elif task_type == "MONITOR":
                self.execute_monitor_task(
                    task, executor, args, kwargs, var_args, var_kwargs
                )
            elif task_type == "NORMAL":
                self.execute_normal_task(
                    task,
                    executor,
                    args,
                    kwargs,
                    var_args,
                    var_kwargs,
                    continue_workgraph,
                )
            else:
                self.process.report(f"Unknown task type {task_type}")
                return self.process.exit_codes.UNKNOWN_TASK_TYPE

    def execute_node_task(
        self, name, executor, kwargs, var_args, var_kwargs, continue_workgraph
    ):
        """Execute a NODE task."""
        results = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
        self.set_task_state_info(name, "process", results)
        self.update_task_state(name)
        if continue_workgraph:
            self.continue_workgraph()

    def execute_data_task(self, name, executor, args, kwargs, continue_workgraph):
        """Execute a DATA task."""
        from aiida_workgraph.utils import create_data_node

        for key in self.ctx._tasks[name]["args"]:
            kwargs.pop(key, None)
        results = create_data_node(executor, args, kwargs)
        self.set_task_state_info(name, "process", results)
        self.update_task_state(name)
        self.ctx._new_data[name] = results
        if continue_workgraph:
            self.continue_workgraph()

    def execute_function_task(
        self, name, executor, kwargs, var_kwargs, continue_workgraph
    ):
        """Execute a CalcFunction or WorkFunction task."""
        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": name})
        try:
            # since aiida 2.5.0, we need to use args_dict to pass the args to the run_get_node
            if var_kwargs is None:
                results, process = run_get_node(executor, **kwargs)
            else:
                results, process = run_get_node(executor, **kwargs, **var_kwargs)
            process.label = name
            self.set_task_state_info(name, "process", process)
            self.update_task_state(name)
        except Exception as e:
            self.logger.error(f"Error in task {name}: {e}")
            self.update_task_state(name, success=False)
        # exclude the current tasks from the next run
        if continue_workgraph:
            self.continue_workgraph()

    def execute_process_task(self, name, executor, kwargs):
        """Execute a CalcJob or WorkChain task."""
        # process = run_get_node(executor, *args, **kwargs)
        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": name})
        # transfer the args to kwargs
        if self.get_task_state_info(name, "action").upper() == "PAUSE":
            self.set_task_state_info(name, "action", "")
            self.process.report(f"Task {name} is created and paused.")
            process = create_and_pause_process(
                self.runner,
                executor,
                kwargs,
                state_msg="Paused through WorkGraph",
            )
            self.set_task_state_info(name, "state", "CREATED")
            process = process.node
        else:
            process = self.process.submit(executor, **kwargs)
            self.set_task_state_info(name, "state", "RUNNING")
        process.label = name
        self.set_task_state_info(name, "process", process)
        self.awaitable_manager.to_context(**{name: process})

    def execute_graph_builder_task(self, task, executor, kwargs, var_args, var_kwargs):
        """Execute a GraphBuilder task."""
        name = task["name"]
        wg = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
        wg.name = name
        wg.group_outputs = self.ctx._tasks[name]["metadata"]["group_outputs"]
        wg.parent_uuid = self.process.node.uuid
        inputs = wg.prepare_inputs(metadata={"call_link_label": name})
        process = self.process.submit(self.process.__class__, inputs=inputs)
        self.set_task_state_info(name, "process", process)
        self.set_task_state_info(name, "state", "RUNNING")
        self.awaitable_manager.to_context(**{name: process})

    def execute_workgraph_task(self, task, kwargs):
        from .utils import prepare_for_workgraph_task

        name = task["name"]
        inputs, _ = prepare_for_workgraph_task(task, kwargs)
        process = self.process.submit(self.process.__class__, inputs=inputs)
        self.set_task_state_info(name, "process", process)
        self.set_task_state_info(name, "state", "RUNNING")
        self.awaitable_manager.to_context(**{name: process})

    def execute_python_job_task(self, task, kwargs, var_kwargs):
        """Execute a PythonJob task."""
        from aiida_pythonjob import PythonJob
        from .utils import prepare_for_python_task

        name = task["name"]
        inputs = prepare_for_python_task(task, kwargs, var_kwargs)
        # since aiida 2.5.0, we can pass inputs directly to the submit, no need to use **inputs
        if self.get_task_state_info(name, "action").upper() == "PAUSE":
            self.set_task_state_info(name, "action", "")
            self.process.report(f"Task {name} is created and paused.")
            process = create_and_pause_process(
                self.runner,
                PythonJob,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            self.set_task_state_info(name, "state", "CREATED")
            process = process.node
        else:
            process = self.process.submit(PythonJob, **inputs)
            self.set_task_state_info(name, "state", "RUNNING")
        process.label = name
        self.set_task_state_info(name, "process", process)
        self.awaitable_manager.to_context(**{name: process})

    def execute_shell_job_task(self, task, kwargs):
        """Execute a ShellJob task."""
        from aiida_shell.calculations.shell import ShellJob
        from .utils import prepare_for_shell_task

        name = task["name"]
        inputs = prepare_for_shell_task(task, kwargs)
        if self.get_task_state_info(name, "action").upper() == "PAUSE":
            self.set_task_state_info(name, "action", "")
            self.process.report(f"Task {name} is created and paused.")
            process = create_and_pause_process(
                self.runner,
                ShellJob,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            self.set_task_state_info(name, "state", "CREATED")
            process = process.node
        else:
            process = self.process.submit(ShellJob, **inputs)
            self.set_task_state_info(name, "state", "RUNNING")
        process.label = name
        self.set_task_state_info(name, "process", process)
        self.awaitable_manager.to_context(**{name: process})

    def execute_while_task(self, task):
        """Execute a WHILE task."""
        # TODO refactor this for while, if and zone
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.are_childen_finished(name)[0]:
            self.update_while_task_state(name)
        else:
            # check the conditions of the while task
            should_run = self.should_run_while_task(name)
            if not should_run:
                self.set_task_state_info(name, "state", "FINISHED")
                self.set_tasks_state(self.ctx._tasks[name]["children"], "SKIPPED")
                self.update_parent_task_state(name)
                self.process.report(
                    f"While Task {name}: Condition not fullilled, task finished. Skip all its children."
                )
            else:
                task["execution_count"] += 1
                self.set_task_state_info(name, "state", "RUNNING")
        self.continue_workgraph()

    def execute_if_task(self, task):
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.are_childen_finished(name)[0]:
            self.update_zone_task_state(name)
        else:
            should_run = self.should_run_if_task(name)
            if should_run:
                self.set_task_state_info(name, "state", "RUNNING")
            else:
                self.set_tasks_state(task["children"], "SKIPPED")
                self.update_zone_task_state(name)
        self.continue_workgraph()

    def execute_zone_task(self, task):
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.are_childen_finished(name)[0]:
            self.update_zone_task_state(name)
        else:
            self.set_task_state_info(name, "state", "RUNNING")
        self.continue_workgraph()

    def execute_get_context_task(self, task, kwargs):
        # get the results from the context
        name = task["name"]
        results = {"result": getattr(self.ctx, kwargs["key"])}
        task["results"] = results
        self.set_task_state_info(name, "state", "FINISHED")
        self.update_parent_task_state(name)
        self.continue_workgraph()

    def execute_set_context_task(self, task, kwargs):
        name = task["name"]
        # get the results from the context
        setattr(self.ctx, kwargs["key"], kwargs["value"])
        self.set_task_state_info(name, "state", "FINISHED")
        self.update_parent_task_state(name)
        self.continue_workgraph()

    def execute_awaitable_task(
        self, task, executor, args, kwargs, var_args, var_kwargs
    ):
        name = task["name"]
        for key in task["args"]:
            kwargs.pop(key, None)
        awaitable_target = asyncio.ensure_future(
            self.run_executor(executor, args, kwargs, var_args, var_kwargs),
            loop=self.process.loop,
        )
        awaitable = self.awaitable_manager.construct_awaitable_function(
            name, awaitable_target
        )
        self.set_task_state_info(name, "state", "RUNNING")
        self.awaitable_manager.to_context(**{name: awaitable})

    def execute_monitor_task(self, task, executor, args, kwargs, var_args, var_kwargs):
        name = task["name"]
        for key in self.ctx._tasks[name]["args"]:
            kwargs.pop(key, None)
        # add function and interval to the args
        args = [
            executor,
            kwargs.pop("interval", 1),
            kwargs.pop("timeout", 3600),
            *args,
        ]
        awaitable_target = asyncio.ensure_future(
            self.run_executor(monitor, args, kwargs, var_args, var_kwargs),
            loop=self.process.loop,
        )
        awaitable = self.awaitable_manager.construct_awaitable_function(
            name, awaitable_target
        )
        self.set_task_state_info(name, "state", "RUNNING")
        # save the awaitable to the temp, so that we can kill it if needed
        self.awaitable_manager.not_persisted_awaitables[name] = awaitable_target
        self.awaitable_manager.to_context(**{name: awaitable})

    def execute_normal_task(
        self, task, executor, args, kwargs, var_args, var_kwargs, continue_workgraph
    ):
        # Normal task is created by decoratoring a function with @task()
        name = task["name"]
        if "context" in task["kwargs"]:
            self.ctx.task_name = name
            kwargs.update({"context": self.ctx})
        for key in self.ctx._tasks[name]["args"]:
            kwargs.pop(key, None)
        try:
            results = self.run_executor(executor, args, kwargs, var_args, var_kwargs)
            self.update_normal_task_state(name, results)
        except Exception as e:
            self.logger.error(f"Error in task {name}: {e}")
            self.update_normal_task_state(name, results=None, success=False)
        if continue_workgraph:
            self.continue_workgraph()

    def get_inputs(
        self, name: str
    ) -> Tuple[
        List[Any],
        Dict[str, Any],
        Optional[List[Any]],
        Optional[Dict[str, Any]],
        Dict[str, Any],
    ]:
        """Get input based on the links."""

        args = []
        args_dict = {}
        kwargs = {}
        var_args = None
        var_kwargs = None
        task = self.ctx._tasks[name]
        properties = task.get("properties", {})
        inputs = {}
        for name, input in task["inputs"].items():
            # print(f"input: {input['name']}")
            if len(input["links"]) == 0:
                inputs[name] = self.ctx_manager.update_context_variable(
                    input["property"]["value"]
                )
            elif len(input["links"]) == 1:
                link = input["links"][0]
                if self.ctx._tasks[link["from_node"]]["results"] is None:
                    inputs[name] = None
                else:
                    # handle the special socket _wait, _outputs
                    if link["from_socket"] == "_wait":
                        continue
                    elif link["from_socket"] == "_outputs":
                        inputs[name] = self.ctx._tasks[link["from_node"]]["results"]
                    else:
                        inputs[name] = get_nested_dict(
                            self.ctx._tasks[link["from_node"]]["results"],
                            link["from_socket"],
                        )
            # handle the case of multiple outputs
            elif len(input["links"]) > 1:
                value = {}
                for link in input["links"]:
                    item_name = f'{link["from_node"]}_{link["from_socket"]}'
                    # handle the special socket _wait, _outputs
                    if link["from_socket"] == "_wait":
                        continue
                    if self.ctx._tasks[link["from_node"]]["results"] is None:
                        value[item_name] = None
                    else:
                        value[item_name] = self.ctx._tasks[link["from_node"]][
                            "results"
                        ][link["from_socket"]]
                inputs[name] = value
        for name in task.get("args", []):
            if name in inputs:
                args.append(inputs[name])
                args_dict[name] = inputs[name]
            else:
                value = self.ctx_manager.update_context_variable(
                    properties[name]["value"]
                )
                args.append(value)
                args_dict[name] = value
        for name in task.get("kwargs", []):
            if name in inputs:
                kwargs[name] = inputs[name]
            else:
                value = self.ctx_manager.update_context_variable(
                    properties[name]["value"]
                )
                kwargs[name] = value
        if task["var_args"] is not None:
            name = task["var_args"]
            if name in inputs:
                var_args = inputs[name]
            else:
                value = self.ctx_manager.update_context_variable(
                    properties[name]["value"]
                )
                var_args = value
        if task["var_kwargs"] is not None:
            name = task["var_kwargs"]
            if name in inputs:
                var_kwargs = inputs[name]
            else:
                value = self.ctx_manager.update_context_variable(
                    properties[name]["value"]
                )
                var_kwargs = value
        return args, kwargs, var_args, var_kwargs, args_dict

    def update_task_state(self, name: str, success=True) -> None:
        """Update task state when the task is finished."""
        task = self.ctx._tasks[name]
        if success:
            node = self.get_task_state_info(name, "process")
            if isinstance(node, ProcessNode):
                # print(f"set task result: {name} process")
                state = node.process_state.value.upper()
                if node.is_finished_ok:
                    self.set_task_state_info(task["name"], "state", state)
                    if task["metadata"]["node_type"].upper() == "WORKGRAPH":
                        # expose the outputs of all the tasks in the workgraph
                        task["results"] = {}
                        outgoing = node.base.links.get_outgoing()
                        for link in outgoing.all():
                            if isinstance(link.node, ProcessNode) and getattr(
                                link.node, "process_state", False
                            ):
                                task["results"][link.link_label] = link.node.outputs
                    else:
                        task["results"] = node.outputs
                        # self.ctx._new_data[name] = task["results"]
                    self.set_task_state_info(task["name"], "state", "FINISHED")
                    self.task_set_context(name)
                    self.process.report(f"Task: {name} finished.")
                # all other states are considered as failed
                else:
                    task["results"] = node.outputs
                    self.on_task_failed(name)
            elif isinstance(node, Data):
                #
                output_name = [
                    output_name
                    for output_name in list(task["outputs"].keys())
                    if output_name not in ["_wait", "_outputs"]
                ][0]
                task["results"] = {output_name: node}
                self.set_task_state_info(task["name"], "state", "FINISHED")
                self.task_set_context(name)
                self.process.report(f"Task: {name} finished.")
            else:
                task.setdefault("results", None)
        else:
            self.on_task_failed(name)
        self.update_parent_task_state(name)

    def on_task_failed(self, name: str) -> None:
        """Handle the case where a task has failed."""
        self.set_task_state_info(name, "state", "FAILED")
        self.set_tasks_state(self.ctx._connectivity["child_node"][name], "SKIPPED")
        self.process.report(f"Task: {name} failed.")
        self.process.error_handler_manager.run_error_handlers(name)

    def update_task(self, task: Task):
        """Update task in the context.
        This is used in error handlers to update the task parameters."""
        tdata = task.to_dict()
        self.ctx._tasks[task.name]["properties"] = tdata["properties"]
        self.ctx._tasks[task.name]["inputs"] = tdata["inputs"]
        self.reset_task(task.name)

    def update_normal_task_state(self, name, results, success=True):
        """Set the results of a normal task.
        A normal task is created by decorating a function with @task().
        """
        from aiida_workgraph.utils import get_sorted_names

        if success:
            task = self.ctx._tasks[name]
            if isinstance(results, tuple):
                # there are two built-in outputs: _wait and _outputs
                if len(task["outputs"]) - 2 != len(results):
                    self.on_task_failed(name)
                    return self.process.exit_codes.OUTPUS_NOT_MATCH_RESULTS
                output_names = get_sorted_names(task["outputs"])[0:-2]
                for i, output_name in enumerate(output_names):
                    task["results"][output_name] = results[i]
            elif isinstance(results, dict):
                task["results"] = results
            else:
                output_name = [
                    output_name
                    for output_name in list(task["outputs"].keys())
                    if output_name not in ["_wait", "_outputs"]
                ][0]
                task["results"][output_name] = results
            self.task_set_context(name)
            self.set_task_state_info(name, "state", "FINISHED")
            self.process.report(f"Task: {name} finished.")
        else:
            self.on_task_failed(name)
        self.update_parent_task_state(name)

    def update_parent_task_state(self, name: str) -> None:
        """Update parent task state."""
        parent_task = self.ctx._tasks[name]["parent_task"]
        if parent_task[0]:
            task_type = self.ctx._tasks[parent_task[0]]["metadata"]["node_type"].upper()
            if task_type == "WHILE":
                self.update_while_task_state(parent_task[0])
            elif task_type == "IF":
                self.update_zone_task_state(parent_task[0])
            elif task_type == "ZONE":
                self.update_zone_task_state(parent_task[0])

    def update_while_task_state(self, name: str) -> None:
        """Update while task state."""
        finished, _ = self.are_childen_finished(name)

        if finished:
            self.process.report(
                f"Wihle Task {name}: this iteration finished. Try to reset for the next iteration."
            )
            # reset the condition tasks
            for link in self.ctx._tasks[name]["inputs"]["conditions"]["links"]:
                self.reset_task(link["from_node"], recursive=False)
            # reset the task and all its children, so that the task can run again
            # do not reset the execution count
            self.reset_task(name, reset_execution_count=False)

    def update_zone_task_state(self, name: str) -> None:
        """Update zone task state."""
        finished, _ = self.are_childen_finished(name)
        if finished:
            self.set_task_state_info(name, "state", "FINISHED")
            self.process.report(f"Task: {name} finished.")
            self.update_parent_task_state(name)

    def should_run_while_task(self, name: str) -> tuple[bool, Any]:
        """Check if the while task should run."""
        # check the conditions of the while task
        not_excess_max_iterations = (
            self.ctx._tasks[name]["execution_count"]
            < self.ctx._tasks[name]["inputs"]["max_iterations"]["property"]["value"]
        )
        conditions = [not_excess_max_iterations]
        _, kwargs, _, _, _ = self.get_inputs(name)
        if isinstance(kwargs["conditions"], list):
            for condition in kwargs["conditions"]:
                value = get_nested_dict(self.ctx, condition)
                conditions.append(value)
        elif isinstance(kwargs["conditions"], dict):
            for _, value in kwargs["conditions"].items():
                conditions.append(value)
        else:
            conditions.append(kwargs["conditions"])
        return False not in conditions

    def should_run_if_task(self, name: str) -> tuple[bool, Any]:
        """Check if the IF task should run."""
        _, kwargs, _, _, _ = self.get_inputs(name)
        flag = kwargs["conditions"]
        if kwargs["invert_condition"]:
            return not flag
        return flag

    def are_childen_finished(self, name: str) -> tuple[bool, Any]:
        """Check if the child tasks are finished."""
        task = self.ctx._tasks[name]
        finished = True
        for name in task["children"]:
            if self.get_task_state_info(name, "state") not in [
                "FINISHED",
                "SKIPPED",
                "FAILED",
            ]:
                finished = False
                break
        return finished, None

    def run_executor(
        self,
        executor: Callable,
        args: List[Any],
        kwargs: Dict[str, Any],
        var_args: Optional[List[Any]],
        var_kwargs: Optional[Dict[str, Any]],
    ) -> Any:
        if var_kwargs is None:
            return executor(*args, **kwargs)
        else:
            return executor(*args, **kwargs, **var_kwargs)

    def apply_task_actions(self, msg: dict) -> None:
        """Apply task actions to the workgraph."""
        action = msg["action"]
        tasks = msg["tasks"]
        self.process.report(f"Action: {action}. {tasks}")
        if action.upper() == "RESET":
            for name in tasks:
                self.reset_task(name)
        elif action.upper() == "PAUSE":
            for name in tasks:
                self.pause_task(name)
        elif action.upper() == "PLAY":
            for name in tasks:
                self.play_task(name)
        elif action.upper() == "SKIP":
            for name in tasks:
                self.skip_task(name)
        elif action.upper() == "KILL":
            for name in tasks:
                self.kill_task(name)

    def pause_task(self, name: str) -> None:
        """Pause task."""
        self.set_task_state_info(name, "action", "PAUSE")
        self.process.report(f"Task {name} action: PAUSE.")

    def play_task(self, name: str) -> None:
        """Play task."""
        self.set_task_state_info(name, "action", "")
        self.process.report(f"Task {name} action: PLAY.")

    def skip_task(self, name: str) -> None:
        """Skip task."""
        self.set_task_state_info(name, "state", "SKIPPED")
        self.process.report(f"Task {name} action: SKIP.")

    def kill_task(self, name: str) -> None:
        """Kill task.
        This is used to kill the awaitable and monitor task.
        """
        if self.get_task_state_info(name, "state") in ["RUNNING"]:
            if self.ctx._tasks[name]["metadata"]["node_type"].upper() in [
                "AWAITABLE",
                "MONITOR",
            ]:
                try:
                    self.awaitable_manager.not_persisted_awaitables[name].cancel()
                    self.set_task_state_info(name, "state", "KILLED")
                    self.process.report(f"Task {name} action: KILLED.")
                except Exception as e:
                    self.logger.error(f"Error in killing task {name}: {e}")

    def check_while_conditions(self) -> bool:
        """Check while conditions.
        Run all condition tasks and check if all the conditions are True.
        """
        self.process.report("Check while conditions.")
        if self.ctx._execution_count >= self.ctx._max_iteration:
            self.process.report("Max iteration reached.")
            return False
        condition_tasks = []
        for c in self.ctx._workgraph["conditions"]:
            task_name, socket_name = c.split(".")
            if "task_name" != "context":
                condition_tasks.append(task_name)
                self.reset_task(task_name)
        self.run_tasks(condition_tasks, continue_workgraph=False)
        conditions = []
        for c in self.ctx._workgraph["conditions"]:
            task_name, socket_name = c.split(".")
            if task_name == "context":
                conditions.append(self.ctx[socket_name])
            else:
                conditions.append(self.ctx._tasks[task_name]["results"][socket_name])
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_tasks_state(condition_tasks, "SKIPPED")
        return should_run

    def check_for_conditions(self) -> bool:
        condition_tasks = [c[0] for c in self.ctx._workgraph["conditions"]]
        self.run_tasks(condition_tasks)
        conditions = [self.ctx._count < len(self.ctx._sequence)] + [
            self.ctx._tasks[c[0]]["results"][c[1]]
            for c in self.ctx._workgraph["conditions"]
        ]
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_tasks_state(condition_tasks, "SKIPPED")
            self.ctx["i"] = self.ctx._sequence[self.ctx._count]
        self.ctx._count += 1
        return should_run

    def reset(self) -> None:
        self.ctx._execution_count += 1
        self.set_tasks_state(self.ctx._tasks.keys(), "PLANNED")
        self.ctx._executed_tasks = []
