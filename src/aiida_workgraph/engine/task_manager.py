from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Callable
from aiida_workgraph.task import Task
from aiida_workgraph.utils import get_nested_dict
import asyncio
from aiida.engine.processes.exit_code import ExitCode
from aiida_workgraph.executors.monitors import monitor
from .task_state import TaskStateManager
from .task_actions import TaskActionManager
import traceback

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
        """
        :param ctx_manager: The object managing the 'ctx' dictionary.
        :param logger: A logger instance.
        :param runner: An AiiDA runner.
        :param process: The AiiDA process object that orchestrates the entire WorkGraph.
        :param awaitable_manager: Manages the global context variables.
        """
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger
        self.runner = runner
        self.process = process
        self.awaitable_manager = awaitable_manager
        # Sub-managers
        self.state_manager = TaskStateManager(
            ctx_manager, logger, process, awaitable_manager
        )
        self.action_manager = TaskActionManager(self.state_manager, logger, process)

    def get_task_executor_dict(self, name: str):
        # if task is a mappped task, we need to get the executor from the parent task
        # TODO: recursive get the parent task
        if "map_data" in self.ctx._tasks[name]:
            parent_task_name = self.ctx._tasks[name]["map_data"]["parent"]
            executor_dict = self.process.node.task_executors[parent_task_name]
        else:
            executor_dict = self.process.node.task_executors[name]
        return executor_dict

    def get_task(self, name: str):
        """Get task from the context."""
        executor_dict = self.get_task_executor_dict(name)
        tdata = {"executor": executor_dict}
        tdata.update(self.ctx._tasks[name])
        task = Task.from_dict(tdata)
        action = self.state_manager.get_task_runtime_info(name, "action").upper()
        task.action = action
        # update task results
        # namespace socket does not have a value, but _value
        for output in task.outputs:
            output.value = get_nested_dict(
                self.ctx._tasks[name]["results"],
                output._name,
                default=output.value if hasattr(output, "value") else None,
            )
        return task

    def set_task_results(self) -> None:
        for name, task in self.ctx._tasks.items():
            if (
                self.state_manager.get_task_runtime_info(name, "action").upper()
                == "RESET"
            ):
                self.state_manager.reset_task(task["name"])
            self.state_manager.update_task_state(name)

    def is_workgraph_finished(self) -> bool:
        """Check if the workgraph is finished.
        For `while` workgraph, we need check its conditions"""
        is_finished = True
        failed_tasks = []
        not_finished_tasks = []
        for name, task in self.ctx._tasks.items():
            # if the task is in mapped state, we need to check its children (mapped tasks)
            if self.state_manager.get_task_runtime_info(task["name"], "state") in [
                "MAPPED"
            ]:
                self.state_manager.update_template_task_state(name)
            elif self.state_manager.get_task_runtime_info(task["name"], "state") in [
                "RUNNING",
                "CREATED",
                "PLANNED",
                "READY",
            ]:
                not_finished_tasks.append(name)
                is_finished = False
            elif (
                self.state_manager.get_task_runtime_info(task["name"], "state")
                == "FAILED"
            ):
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
        """
        Resume the WorkGraph by looking for tasks that are ready to run.
        """
        self.process.report("Continue workgraph.")
        task_to_run = []
        for name, task in self.ctx._tasks.items():
            # update task state
            if (
                self.state_manager.get_task_runtime_info(task["name"], "state")
                in [
                    "CREATED",
                    "RUNNING",
                    "FINISHED",
                    "FAILED",
                    "SKIPPED",
                    "MAPPED",
                ]
                or name in self.ctx._executed_tasks
            ):
                continue
            ready, _ = self.state_manager.is_task_ready_to_run(name)
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
        from aiida_workgraph.utils import update_nested_dict_with_special_keys
        from node_graph.executor import NodeExecutor

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
            if (
                name in self.ctx._executed_tasks
                or self.state_manager.get_task_runtime_info(name, "state")
                in ["SKIPPED"]
            ):
                continue
            self.ctx._executed_tasks.append(name)
            print("-" * 60)

            self.process.report(
                f"Run task: {name}, type: {task['metadata']['node_type']}"
            )
            executor_dict = self.get_task_executor_dict(name)
            executor = NodeExecutor(**executor_dict).executor if executor_dict else None
            args, kwargs, var_args, var_kwargs, args_dict = self.get_inputs(name)
            for i, key in enumerate(self.ctx._tasks[name]["args"]):
                kwargs[key] = args[i]
            # update the port namespace
            kwargs = update_nested_dict_with_special_keys(kwargs)

            print("kwargs: ", kwargs)
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
                self.execute_function_task(name, kwargs, var_kwargs, continue_workgraph)
            elif task_type in [
                "CALCJOB",
                "WORKCHAIN",
                "SHELLJOB",
                "PYTHONJOB",
                "WORKGRAPH",
                "GRAPH_BUILDER",
            ]:
                self.execute_process_task(
                    name, args=args, kwargs=kwargs, var_kwargs=var_kwargs
                )
            elif task_type == "WHILE":
                self.execute_while_task(task)
            elif task_type == "IF":
                self.execute_if_task(task)
            elif task_type == "ZONE":
                self.execute_zone_task(task)
            elif task_type == "MAP":
                self.execute_map_task(task, kwargs)
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
        self.state_manager.set_task_runtime_info(name, "process", results)
        self.state_manager.update_task_state(name)
        if continue_workgraph:
            self.continue_workgraph()

    def execute_data_task(self, name, executor, args, kwargs, continue_workgraph):
        """Execute a DATA task."""
        from aiida_workgraph.utils import create_data_node

        for key in self.ctx._tasks[name]["args"]:
            kwargs.pop(key, None)
        results = create_data_node(executor, args, kwargs)
        self.state_manager.set_task_runtime_info(name, "process", results)
        self.state_manager.update_task_state(name)
        self.ctx._new_data[name] = results
        if continue_workgraph:
            self.continue_workgraph()

    def execute_function_task(self, name, kwargs, var_kwargs, continue_workgraph):
        """Execute a CalcFunction or WorkFunction task."""

        try:
            task = self.get_task(name)
            process, _ = task.execute(None, kwargs, var_kwargs)
            self.state_manager.set_task_runtime_info(name, "process", process)
            self.state_manager.update_task_state(name)
        except Exception as e:
            error_traceback = traceback.format_exc()  # Capture the full traceback
            self.logger.error(
                f"Error in task {name}: {e}\n{error_traceback}"
            )  # Log the error with traceback
            self.state_manager.update_task_state(name, success=False)
        # exclude the current tasks from the next run
        if continue_workgraph:
            self.continue_workgraph()

    def execute_process_task(self, name, args=None, kwargs=None, var_kwargs=None):
        """Execute a CalcJob or WorkChain task."""
        try:
            task = self.get_task(name)
            process, state = task.execute(
                engine_process=self.process,
                args=args,
                kwargs=kwargs,
                var_kwargs=var_kwargs,
            )
            self.state_manager.set_task_runtime_info(name, "state", state)
            self.state_manager.set_task_runtime_info(name, "action", "")
            self.state_manager.set_task_runtime_info(name, "process", process)
            # update the parent task state of mappped tasks
            if "map_data" in self.ctx._tasks[name]:
                parent_task_name = self.ctx._tasks[name]["map_data"]["parent"]
                if self.process.node.get_task_state(parent_task_name) == "PLANNED":
                    self.process.node.set_task_state(parent_task_name, state)
            self.awaitable_manager.to_context(**{name: process})
        except Exception as e:
            error_traceback = traceback.format_exc()  # Capture the full traceback
            self.logger.error(
                f"Error in task {name}: {e}\n{error_traceback}"
            )  # Log the error with traceback
            self.state_manager.update_task_state(name, success=False)

    def execute_while_task(self, task):
        """Execute a WHILE task."""
        # TODO refactor this for while, if and zone
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.state_manager.are_childen_finished(name)[0]:
            self.state_manager.update_while_task_state(name)
        else:
            # check the conditions of the while task
            should_run = self.should_run_while_task(name)
            if not should_run:
                self.state_manager.set_task_runtime_info(name, "state", "FINISHED")
                self.state_manager.set_tasks_state(
                    self.ctx._tasks[name]["children"], "SKIPPED"
                )
                self.state_manager.update_parent_task_state(name)
                self.process.report(
                    f"While Task {name}: Condition not fullilled, task finished. Skip all its children."
                )
            else:
                task["execution_count"] += 1
                self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
        self.continue_workgraph()

    def execute_if_task(self, task):
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.state_manager.are_childen_finished(name)[0]:
            self.state_manager.update_zone_task_state(name)
        else:
            should_run = self.should_run_if_task(name)
            if should_run:
                self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
            else:
                self.state_manager.set_tasks_state(task["children"], "SKIPPED")
                self.state_manager.update_zone_task_state(name)
        self.continue_workgraph()

    def execute_zone_task(self, task):
        # in case of an empty zone, it will finish immediately
        name = task["name"]
        if self.state_manager.are_childen_finished(name)[0]:
            self.state_manager.update_zone_task_state(name)
        else:
            self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
        self.continue_workgraph()

    def execute_map_task(self, task, kwargs):
        """
        1. Clone the subgraph tasks for each item in `source`.
        2. Mark this MAP node as running and schedule a continuation.
        """
        name = task["name"]
        if self.state_manager.are_childen_finished(name)[0]:
            self.state_manager.update_zone_task_state(name)
        else:
            self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
            source = kwargs["source"]
            placeholder = kwargs.get("placeholder", "map_input")
            for prefix, value in source.items():
                self.generate_mapped_tasks(
                    task, prefix=prefix, placeholder=placeholder, value=value
                )

        self.continue_workgraph()

    def execute_get_context_task(self, task, kwargs):
        # get the results from the context
        name = task["name"]
        results = {"result": getattr(self.ctx, kwargs["key"])}
        task["results"] = results
        self.state_manager.set_task_runtime_info(name, "state", "FINISHED")
        self.state_manager.update_parent_task_state(name)
        self.continue_workgraph()

    def execute_set_context_task(self, task, kwargs):
        name = task["name"]
        # get the results from the context
        setattr(self.ctx, kwargs["key"], kwargs["value"])
        self.state_manager.set_task_runtime_info(name, "state", "FINISHED")
        self.state_manager.update_parent_task_state(name)
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
        self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
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
        self.state_manager.set_task_runtime_info(name, "state", "RUNNING")
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
            self.state_manager.update_normal_task_state(name, results)
        except Exception as e:
            self.logger.error(f"Error in task {name}: {e}")
            self.state_manager.update_normal_task_state(
                name, results=None, success=False
            )
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
        from node_graph.utils import collect_values_inside_namespace

        args = []
        args_dict = {}
        kwargs = {}
        var_args = None
        var_kwargs = None
        task = self.ctx._tasks[name]
        inputs = {}
        for name, prop in task.get("properties", {}).items():
            inputs[name] = self.ctx_manager.update_context_variable(prop["value"])
        for name, input in task["inputs"].items():
            if input["identifier"] == "workgraph.namespace":
                # inputs[name] = self.ctx_manager.update_context_variable(input["value"])
                value = collect_values_inside_namespace(input)
                if value:
                    inputs[name] = value
            else:
                value = self.ctx_manager.update_context_variable(
                    input["property"]["value"]
                )
                if value is not None:
                    inputs[name] = value
        for name, links in task["input_links"].items():
            if len(links) == 1:
                link = links[0]
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
            elif len(links) > 1:
                value = {}
                for link in links:
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
        for name, input in inputs.items():
            # only need to check the top level key
            key = name.split(".")[0]
            if key in task["args"]:
                args.append(input)
                args_dict[name] = input
            elif key in task["kwargs"]:
                kwargs[name] = input
            elif key == task["var_args"]:
                var_args = input
            elif key == task["var_kwargs"]:
                var_kwargs = input
        return args, kwargs, var_args, var_kwargs, args_dict

    def update_task(self, task: Task):
        """Update task in the context.
        This is used in error handlers to update the task parameters."""
        tdata = task.to_dict()
        self.ctx._tasks[task.name]["properties"] = tdata["properties"]
        self.ctx._tasks[task.name]["inputs"] = tdata["inputs"]
        self.state_manager.reset_task(task.name)

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
                self.state_manager.reset_task(task_name)
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
            self.state_manager.set_tasks_state(condition_tasks, "SKIPPED")
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
            self.state_manager.set_tasks_state(condition_tasks, "SKIPPED")
            self.ctx["i"] = self.ctx._sequence[self.ctx._count]
        self.ctx._count += 1
        return should_run

    def reset(self) -> None:
        self.ctx._execution_count += 1
        self.state_manager.set_tasks_state(self.ctx._tasks.keys(), "PLANNED")
        self.ctx._executed_tasks = []

    def get_all_children(self, name: str) -> List[str]:
        """Find all children of the zone_task, and their children recursively"""
        child_tasks = {}
        task = self.ctx._tasks[name]
        for child_task in task.get("children", []):
            child_tasks[child_task] = self.get_all_children(child_task)
            children = self.ctx._tasks[child_task].get("children")
            if children:
                child_tasks.update(self.get_all_children(child_task))
        return child_tasks

    def generate_mapped_tasks(
        self, zone_task: dict, prefix: str, value: any, placeholder: str = "map_input"
    ) -> None:
        """
        Recursively clone the subgraph starting from zone_children,
        rewriting references to old tasks with new task names.
        """
        import uuid

        new_tasks = {}
        child_tasks = self.get_all_children(zone_task["name"])
        for child_task in child_tasks:
            # since the child task is mapped, it should be skipped
            self.state_manager.set_task_runtime_info(child_task, "state", "MAPPED")
            # keep track of the mapped tasks
            self.ctx._tasks[child_task].setdefault("mapped_tasks", {})
            new_name = f"{prefix}_{child_task}"
            new_task_data = self.copy_task(self.ctx._tasks[child_task])
            new_task_data["name"] = new_name
            new_task_data["map_data"] = {"parent": child_task, "prefix": prefix}
            new_task_data["uuid"] = str(uuid.uuid4())
            # Reset runtime states
            new_task_data["results"] = {}
            self.state_manager.set_task_runtime_info(new_name, "state", "PLANNED")
            self.state_manager.set_task_runtime_info(new_name, "action", "")
            # Insert new_data in ctx._tasks
            self.ctx._tasks[new_name] = new_task_data
            self.ctx._tasks[child_task]["mapped_tasks"][prefix] = new_task_data
            new_tasks[child_task] = new_task_data

        # fix references in the newly mapped tasks (children, input_links, etc.)
        self._patch_cloned_tasks(new_tasks, value, placeholder=placeholder)

        # update ctx._connectivity so the new tasks are recognized in child_node, zone references, etc.
        self._patch_connectivity(new_tasks)

    def copy_task(self, task: dict) -> dict:
        from aiida_workgraph.utils import shallow_copy_nested_dict
        from copy import deepcopy

        new_task_data = shallow_copy_nested_dict(task)
        new_task_data["input_links"] = deepcopy(task["input_links"])
        return new_task_data

    def _patch_cloned_tasks(
        self,
        new_tasks: dict[str, str],
        value: any,
        placeholder: str = "map_input",
    ):
        """
        For each newly mapped task, fix references (children, input_links, etc.)
        from old_name -> new_name.
        """
        for task in new_tasks.values():
            # update input
            for input in task["inputs"].values():
                if "property" not in input:
                    continue
                if (
                    isinstance(input["property"]["value"], str)
                    and placeholder in input["property"]["value"]
                ):
                    input["property"]["value"] = value
            # fix children references
            new_children = []
            for child in task.get("children", []):
                new_children.append(new_tasks[child]["name"])
            task["children"] = new_children
            # since this is a newly created task, it should not have any mapped tasks
            task.pop("mapped_tasks", None)
            # fix parent reference
            new_parent = []
            for parent in task.get("parent_task", []):
                if parent in new_tasks:
                    new_parent.append(new_tasks[parent]["name"])
                else:
                    new_parent.append(parent)
            task["parent_task"] = new_parent
            # fix input_links references
            for links in task["input_links"].values():
                for link in links:
                    # this node is mapped, so we need to update the link
                    if link["from_node"] in new_tasks:
                        link["from_node"] = new_tasks[link["from_node"]]["name"]

    def _patch_connectivity(self, new_tasks: dict[str, str]):
        """
        Update the global connectivity for newly created tasks.
        """
        for name, task in new_tasks.items():
            # child_node
            new_child_node = []
            for child_task in self.ctx._connectivity["child_node"][name]:
                if child_task in new_tasks:
                    new_child_node.append(new_tasks[child_task]["name"])
                else:
                    new_child_node.append(child_task)
            self.ctx._connectivity["child_node"][task["name"]] = new_child_node
            # input_tasks
            new_input_tasks = []
            for input_task in self.ctx._connectivity["zone"][name]["input_tasks"]:
                if input_task in new_tasks:
                    new_input_tasks.append(new_tasks[input_task]["name"])
                else:
                    new_input_tasks.append(input_task)
            self.ctx._connectivity["zone"][task["name"]] = {
                "input_tasks": new_input_tasks
            }
