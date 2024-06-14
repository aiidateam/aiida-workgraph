"""AiiDA workflow components: WorkGraph."""
from __future__ import annotations

import collections.abc
import functools
import logging
import typing as t

from plumpy.persistence import auto_persist
from plumpy.process_states import Continue, Wait
from plumpy.workchains import _PropagateReturn

from aiida.common import exceptions
from aiida.common.extendeddicts import AttributeDict
from aiida.common.lang import override
from aiida import orm
from aiida.orm import load_node, Node, ProcessNode, WorkChainNode

from aiida.engine.processes.exit_code import ExitCode
from aiida.engine.processes.process import Process

from aiida.engine.processes.workchains.awaitable import (
    Awaitable,
    AwaitableAction,
    AwaitableTarget,
    construct_awaitable,
)
from aiida.engine.processes.workchains.workchain import Protect, WorkChainSpec
from aiida.engine import run_get_node


if t.TYPE_CHECKING:
    from aiida.engine.runners import Runner  # pylint: disable=unused-import

__all__ = "WorkGraph"


MAX_NUMBER_AWAITABLES_MSG = "The maximum number of subprocesses has been reached: {}. Cannot launch the job: {}."


@auto_persist("_awaitables")
class WorkGraphEngine(Process, metaclass=Protect):
    """The `WorkGraph` class is used to construct workflows in AiiDA."""

    # used to create a process node that represents what happened in this process.
    _node_class = WorkChainNode
    _spec_class = WorkChainSpec
    _CONTEXT = "CONTEXT"

    def __init__(
        self,
        inputs: dict | None = None,
        logger: logging.Logger | None = None,
        runner: "Runner" | None = None,
        enable_persistence: bool = True,
    ) -> None:
        """Construct a WorkGraph instance.

        :param inputs: work graph inputs
        :param logger: aiida logger
        :param runner: work graph runner
        :param enable_persistence: whether to persist this work graph

        """

        super().__init__(inputs, logger, runner, enable_persistence=enable_persistence)

        self._awaitables: list[Awaitable] = []
        self._context = AttributeDict()

    @classmethod
    def define(cls, spec: WorkChainSpec) -> None:
        super().define(spec)
        spec.input("input_file", valid_type=orm.SinglefileData, required=False)
        spec.input_namespace(
            "wg", dynamic=True, required=False, help="WorkGraph inputs"
        )
        spec.input_namespace("input_tasks", dynamic=True, required=False)
        spec.exit_code(2, "ERROR_SUBPROCESS", message="A subprocess has failed.")

        spec.output_namespace("new_data", dynamic=True)
        spec.output_namespace("group_outputs", dynamic=True)
        spec.output(
            "execution_count",
            valid_type=orm.Int,
            required=False,
            help="The number of time the WorkGraph runs.",
        )
        #
        spec.exit_code(
            201, "UNKNOWN_MESSAGE_TYPE", message="The message type is unknown."
        )
        spec.exit_code(202, "UNKNOWN_TASK_TYPE", message="The task type is unknown.")
        #
        spec.exit_code(
            301,
            "OUTPUS_NOT_MATCH_RESULTS",
            message="The outputs of the process do not match the results.",
        )
        spec.exit_code(
            302,
            "TASK_FAILED",
            message="Some of the tasks failed.",
        )
        spec.exit_code(
            303,
            "TASK_NON_ZERO_EXIT_STATUS",
            message="Some of the tasks exited with non-zero status.",
        )

    @property
    def ctx(self) -> AttributeDict:
        """Get the context."""
        return self._context

    @override
    def save_instance_state(
        self, out_state: t.Dict[str, t.Any], save_context: t.Any
    ) -> None:
        """Save instance state.

        :param out_state: state to save in

        :param save_context:
        :type save_context: :class:`!plumpy.persistence.LoadSaveContext`

        """
        super().save_instance_state(out_state, save_context)
        # Save the context
        out_state[self._CONTEXT] = self.ctx

    @override
    def load_instance_state(
        self, saved_state: t.Dict[str, t.Any], load_context: t.Any
    ) -> None:
        super().load_instance_state(saved_state, load_context)
        # Load the context
        self._context = saved_state[self._CONTEXT]

        self.set_logger(self.node.logger)

        if self._awaitables:
            # this is a new runner, so we need to re-register the callbacks
            self.ctx._awaitable_actions = []
            self._action_awaitables()

    def _resolve_nested_context(self, key: str) -> tuple[AttributeDict, str]:
        """
        Returns a reference to a sub-dictionary of the context and the last key,
        after resolving a potentially segmented key where required sub-dictionaries are created as needed.

        :param key: A key into the context, where words before a dot are interpreted as a key for a sub-dictionary
        """
        ctx = self.ctx
        ctx_path = key.split(".")

        for index, path in enumerate(ctx_path[:-1]):
            try:
                ctx = ctx[path]
            except KeyError:  # see below why this is the only exception we have to catch here
                ctx[
                    path
                ] = AttributeDict()  # create the sub-dict and update the context
                ctx = ctx[path]
                continue

            # Notes:
            # * the first ctx (self.ctx) is guaranteed to be an AttributeDict, hence the post-"dereference" checking
            # * the values can be many different things: on insertion they are either AtrributeDict, List or Awaitables
            #   (subclasses of AttributeDict) but after resolution of an Awaitable this will be the value itself
            # * assumption: a resolved value is never a plain AttributeDict, on the other hand if a resolved Awaitable
            #   would be an AttributeDict we can append things to it since the order of tasks is maintained.
            if type(ctx) != AttributeDict:  # pylint: disable=C0123
                raise ValueError(
                    f"Can not update the context for key `{key}`: "
                    f' found instance of `{type(ctx)}` at `{".".join(ctx_path[:index + 1])}`, expected AttributeDict'
                )

        return ctx, ctx_path[-1]

    def _insert_awaitable(self, awaitable: Awaitable) -> None:
        """Insert an awaitable that should be terminated before before continuing to the next step.

        :param awaitable: the thing to await
        """
        ctx, key = self._resolve_nested_context(awaitable.key)

        # Already assign the awaitable itself to the location in the context container where it is supposed to end up
        # once it is resolved. This is especially important for the `APPEND` action, since it needs to maintain the
        # order, but the awaitables will not necessarily be resolved in the order in which they are added. By using the
        # awaitable as a placeholder, in the `_resolve_awaitable`, it can be found and replaced by the resolved value.
        if awaitable.action == AwaitableAction.ASSIGN:
            ctx[key] = awaitable
        elif awaitable.action == AwaitableAction.APPEND:
            ctx.setdefault(key, []).append(awaitable)
        else:
            raise AssertionError(f"Unsupported awaitable action: {awaitable.action}")

        self._awaitables.append(
            awaitable
        )  # add only if everything went ok, otherwise we end up in an inconsistent state
        self._update_process_status()

    def _resolve_awaitable(self, awaitable: Awaitable, value: t.Any) -> None:
        """Resolve an awaitable.

        Precondition: must be an awaitable that was previously inserted.

        :param awaitable: the awaitable to resolve
        """
        ctx, key = self._resolve_nested_context(awaitable.key)

        if awaitable.action == AwaitableAction.ASSIGN:
            ctx[key] = value
        elif awaitable.action == AwaitableAction.APPEND:
            # Find the same awaitable inserted in the context
            container = ctx[key]
            for index, placeholder in enumerate(container):
                if (
                    isinstance(placeholder, Awaitable)
                    and placeholder.pk == awaitable.pk
                ):
                    container[index] = value
                    break
            else:
                raise AssertionError(
                    f"Awaitable `{awaitable.pk} was not in `ctx.{awaitable.key}`"
                )
        else:
            raise AssertionError(f"Unsupported awaitable action: {awaitable.action}")

        awaitable.resolved = True
        self._awaitables.remove(
            awaitable
        )  # remove only if everything went ok, otherwise we may lose track

        if not self.has_terminated():
            # the process may be terminated, for example, if the process was killed or excepted
            # then we should not try to update it
            self._update_process_status()

    @Protect.final
    def to_context(self, **kwargs: Awaitable | ProcessNode) -> None:
        """Add a dictionary of awaitables to the context.

        This is a convenience method that provides syntactic sugar, for a user to add multiple intersteps that will
        assign a certain value to the corresponding key in the context of the work graph.
        """
        for key, value in kwargs.items():
            awaitable = construct_awaitable(value)
            awaitable.key = key
            self._insert_awaitable(awaitable)

    def _update_process_status(self) -> None:
        """Set the process status with a message accounting the current sub processes that we are waiting for."""
        if self._awaitables:
            status = f"Waiting for child processes: {', '.join([str(_.pk) for _ in self._awaitables])}"
            self.node.set_process_status(status)
        else:
            self.node.set_process_status(None)

    @override
    def run(self) -> t.Any:
        self.setup()
        return self._do_step()

    def _do_step(self) -> t.Any:
        """Execute the next step in the workgraph and return the result.

        If any awaitables were created, the process will enter in the Wait state,
        otherwise it will go to Continue.
        """
        # we will not remove the awaitables here,
        # we resume the workgraph in the callback function even
        # there are some awaitables left
        # self._awaitables = []
        result: t.Any = None

        try:
            self.continue_workgraph()
        except _PropagateReturn as exception:
            finished, result = True, exception.exit_code
        else:
            finished, result = self.is_workgraph_finished()

        # If the workgraph is finished or the result is an ExitCode, we exit by returning
        if finished:
            if isinstance(result, ExitCode):
                return result
            else:
                return self.finalize()

        if self._awaitables:
            return Wait(self._do_step, "Waiting before next step")

        return Continue(self._do_step)

    def _store_nodes(self, data: t.Any) -> None:
        """Recurse through a data structure and store any unstored nodes that are found along the way

        :param data: a data structure potentially containing unstored nodes
        """
        if isinstance(data, Node) and not data.is_stored:
            data.store()
        elif isinstance(data, collections.abc.Mapping):
            for _, value in data.items():
                self._store_nodes(value)
        elif isinstance(data, collections.abc.Sequence) and not isinstance(data, str):
            for value in data:
                self._store_nodes(value)

    @override
    @Protect.final
    def on_exiting(self) -> None:
        """Ensure that any unstored nodes in the context are stored, before the state is exited

        After the state is exited the next state will be entered and if persistence is enabled, a checkpoint will
        be saved. If the context contains unstored nodes, the serialization necessary for checkpointing will fail.
        """
        super().on_exiting()
        try:
            self._store_nodes(self.ctx)
        except Exception:  # pylint: disable=broad-except
            # An uncaught exception here will have bizarre and disastrous consequences
            self.logger.exception("exception in _store_nodes called in on_exiting")

    @Protect.final
    def on_wait(self, awaitables: t.Sequence[t.Awaitable]):
        """Entering the WAITING state."""
        super().on_wait(awaitables)
        if self._awaitables:
            self._action_awaitables()
        else:
            self.call_soon(self.resume)

    def _action_awaitables(self) -> None:
        """Handle the awaitables that are currently registered with the work chain.

        Depending on the class type of the awaitable's target a different callback
        function will be bound with the awaitable and the runner will be asked to
        call it when the target is completed
        """
        for awaitable in self._awaitables:
            # if the waitable already has a callback, skip
            if awaitable.pk in self.ctx._awaitable_actions:
                continue
            if awaitable.target == AwaitableTarget.PROCESS:
                callback = functools.partial(
                    self.call_soon, self._on_awaitable_finished, awaitable
                )
                self.runner.call_on_process_finish(awaitable.pk, callback)
                self.ctx._awaitable_actions.append(awaitable.pk)
            else:
                assert f"invalid awaitable target '{awaitable.target}'"

    def _on_awaitable_finished(self, awaitable: Awaitable) -> None:
        """Callback function, for when an awaitable process instance is completed.

        The awaitable will be effectuated on the context of the work chain and removed from the internal list. If all
        awaitables have been dealt with, the work chain process is resumed.

        :param awaitable: an Awaitable instance
        """
        print("on awaitable finished: ", awaitable.key)
        self.logger.info(
            "received callback that awaitable %d has terminated", awaitable.pk
        )

        try:
            node = load_node(awaitable.pk)
        except (exceptions.MultipleObjectsError, exceptions.NotExistent):
            raise ValueError(
                f"provided pk<{awaitable.pk}> could not be resolved to a valid Node instance"
            )

        if awaitable.outputs:
            value = {
                entry.link_label: entry.node for entry in node.base.links.get_outgoing()
            }
        else:
            value = node  # type: ignore

        self._resolve_awaitable(awaitable, value)

        # node finished, update the task state and result
        # udpate the task state
        self.update_task_state(awaitable.key)
        # try to resume the workgraph, if the workgraph is already resumed
        # by other awaitable, this will not work
        try:
            self.resume()
        except Exception as e:
            print(e)

    def _build_process_label(self) -> str:
        """Use the workgraph name as the process label."""
        return f"WorkGraph<{self.inputs.wg['name']}>"

    def setup(self) -> None:
        # track if the awaitable callback is added to the runner
        self.ctx._awaitable_actions = []
        self.ctx.new_data = dict()
        self.ctx.input_tasks = dict()
        # read the latest workgraph data
        wgdata = self.read_wgdata_from_base()
        self.init_ctx(wgdata)
        #
        self.ctx.msgs = []
        self.ctx._execution_count = 0
        # init task results
        self.set_task_results()
        # while workgraph
        if self.ctx.workgraph["workgraph_type"].upper() == "WHILE":
            self.ctx._max_iteration = self.ctx.workgraph.get("max_iteration", 1000)
            should_run = self.check_while_conditions()
            if not should_run:
                self.set_task_state(self.ctx.tasks.keys(), "SKIPPED")
        # for workgraph
        if self.ctx.workgraph["workgraph_type"].upper() == "FOR":
            should_run = self.check_for_conditions()
            if not should_run:
                self.set_task_state(self.ctx.tasks.keys(), "SKIPPED")

    def setup_ctx_workgraph(self, wgdata: t.Dict[str, t.Any]) -> None:
        """setup the workgraph in the context."""
        self.ctx.tasks = wgdata["tasks"]
        self.ctx.links = wgdata["links"]
        self.ctx.connectivity = wgdata["connectivity"]
        self.ctx.ctrl_links = wgdata["ctrl_links"]
        self.ctx.workgraph = wgdata

    def read_wgdata_from_base(self) -> t.Dict[str, t.Any]:
        """Read workgraph data from base.extras."""
        from aiida.orm.utils.serialize import deserialize_unsafe

        wgdata = deserialize_unsafe(self.node.base.extras.get("workgraph"))
        return wgdata

    def update_workgraph_from_base(self) -> None:
        """Update the ctx from base.extras."""
        wgdata = self.read_wgdata_from_base()
        for name, task in wgdata["tasks"].items():
            task["state"] = self.ctx.tasks[name]["state"]
            task["results"] = self.ctx.tasks[name].get("results")
            task["process"] = self.ctx.tasks[name].get("process")
        self.setup_ctx_workgraph(wgdata)

    def init_ctx(self, wgdata: t.Dict[str, t.Any]) -> None:
        """Init the context from the workgraph data."""
        from aiida_workgraph.utils import update_nested_dict

        # set up the context variables
        self.ctx.max_number_awaitables = (
            wgdata["max_number_jobs"] if wgdata["max_number_jobs"] else 1000000
        )
        self.ctx["_count"] = 0
        for key, value in wgdata["context"].items():
            key = key.replace("__", ".")
            update_nested_dict(self.ctx, key, value)
        # set up the workgraph
        self.setup_ctx_workgraph(wgdata)

    def set_task_results(self) -> None:
        for _, task in self.ctx.tasks.items():
            if task.get("process"):
                if isinstance(task["process"], str):
                    task["process"] = orm.load_node(task["process"])
                self.set_task_result(task)
            self.set_task_result(task)

    def set_task_result(self, task: t.Dict[str, t.Any]) -> None:
        name = task["name"]
        # print(f"set task result: {name}")
        if task.get("process"):
            # print(f"set task result: {name} process")
            state = task["process"].process_state.value.upper()
            if task["process"].is_finished_ok:
                task["state"] = state
                if task["metadata"]["node_type"].upper() == "GRAPH_BUILDER":
                    # expose the outputs of workgraph
                    task["results"] = getattr(
                        task["process"].outputs, "group_outputs", None
                    )
                    # self.ctx.new_data[name] = outputs
                elif task["metadata"]["node_type"].upper() == "WORKGRAPH":
                    # expose the outputs of all the tasks in the workgraph
                    task["results"] = {}
                    outgoing = task["process"].base.links.get_outgoing()
                    for link in outgoing.all():
                        if isinstance(link.node, ProcessNode) and getattr(
                            link.node, "process_state", False
                        ):
                            task["results"][link.link_label] = link.node.outputs
                else:
                    task["results"] = task["process"].outputs
                    # self.ctx.new_data[name] = task["results"]
                self.ctx.tasks[name]["state"] = "FINISHED"
                self.task_to_context(name)
                self.report(f"Task: {name} finished.")
            # all other states are considered as failed
            else:
                task["state"] = state
                task["results"] = task["process"].outputs
                # self.ctx.new_data[name] = task["results"]
                self.ctx.tasks[name]["state"] = "FAILED"
                # set child tasks state to SKIPPED
                self.set_task_state(
                    self.ctx.connectivity["child_node"][name], "SKIPPED"
                )
                self.report(f"Task: {name} failed.")
        else:
            task["results"] = None

    def apply_actions(self) -> None:
        """Apply actions to the workgraph.
        The actions are stored in the base.extras["workgraph_queue"].
        The index of the last applied action is stored in the base.extras["workgraph_queue_index"].
        """
        msgs = self.node.base.extras.get("workgraph_queue", [])
        index = self.node.base.extras.get("workgraph_queue_index", 0)
        for msg in msgs[index:]:
            header, msg = msg.split(",")
            if header == "task":
                self.apply_task_actions(msg)
            else:
                self.report(f"Unknow message type {msg}")
            index += 1
            self.report("Apply actions: {}".format(msg))
            msgs = self.node.base.extras.set("workgraph_queue_index", index)

    def apply_task_actions(self, msg: str) -> None:
        """Apply task actions to the workgraph."""
        name, action = msg.split(":")
        print("apply task actions: ", name, action)
        if action == "RESET":
            self.reset_task(name)

    def reset_task(self, name: str) -> None:
        """Reset task."""
        self.ctx.tasks[name]["state"] = "CREATED"
        # reset its child tasks
        names = self.ctx.connectivity["child_node"][name]
        for name in names:
            self.ctx.tasks[name]["state"] = "CREATED"
            self.ctx.tasks[name]["result"] = None
            self.ctx.tasks[name]["process"] = None

    def continue_workgraph(self, exclude: t.Optional[t.List[str]] = None) -> None:
        print("Continue workgraph.")
        exclude = exclude or []
        self.report("Continue workgraph.")
        self.update_workgraph_from_base()
        self.apply_actions()
        task_to_run = []
        for name, task in self.ctx.tasks.items():
            # update task state
            if task["state"] in ["RUNNING", "FINISHED", "FAILED", "SKIPPED"]:
                continue
            ready, output = self.check_parent_state(name)
            if ready and name not in exclude:
                task_to_run.append(name)
        #
        self.report("tasks ready to run: {}".format(",".join(task_to_run)))
        self.run_taskss(task_to_run)

    def update_task_state(self, name: str) -> None:
        """Update task state if task is a Awaitable."""
        print("update task state: ", name)
        task = self.ctx.tasks[name]
        if (
            task["metadata"]["node_type"].upper()
            in [
                "CALCFUNCTION",
                "WORKFUNCTION",
                "CALCJOB",
                "WORKCHAIN",
                "GRAPH_BUILDER",
                "WORKGRAPH",
                "PYTHONTASK",
                "SHELLTASK",
            ]
            and task["state"] == "RUNNING"
        ):
            self.set_task_result(task)

    def is_workgraph_finished(self) -> bool:
        """Check if the workgraph is finished.
        For `while` workgraph, we need check its conditions"""
        is_finished = True
        failed_tasks = []
        for name, task in self.ctx.tasks.items():
            # self.update_task_state(name)
            print("task: ", name, task["state"])
            if task["state"] in ["RUNNING", "CREATED", "READY"]:
                is_finished = False
            elif task["state"] == "FAILED":
                failed_tasks.append(name)
        if is_finished:
            if self.ctx.workgraph["workgraph_type"].upper() == "WHILE":
                should_run = self.check_while_conditions()
                is_finished = not should_run
            if self.ctx.workgraph["workgraph_type"].upper() == "FOR":
                should_run = self.check_for_conditions()
                is_finished = not should_run
        print("is workgraph finished: ", is_finished)
        if is_finished and len(failed_tasks) > 0:
            message = f"WorkGraph finished, but tasks: {failed_tasks} failed. Thus all their child tasks are skipped."
            self.report(message)
            result = ExitCode(302, message)
        else:
            result = None
        return is_finished, result

    def check_while_conditions(self) -> bool:
        """Check while conditions.
        Run all condition tasks and check if all the conditions are True.
        """
        print("Is a while workgraph")
        print(
            "execution count: ",
            self.ctx._execution_count,
            "max iteration: ",
            self.ctx._max_iteration,
        )
        self.report("Check while conditions.")
        if self.ctx._execution_count >= self.ctx._max_iteration:
            print("Max iteration reached.")
            self.report("Max iteration reached.")
            return False
        condition_tasks = []
        for c in self.ctx.workgraph["conditions"]:
            task_name, socket_name = c.split(".")
            if "task_name" != "context":
                condition_tasks.append(task_name)
        self.run_taskss(condition_tasks, continue_workgraph=False)
        conditions = []
        for c in self.ctx.workgraph["conditions"]:
            task_name, socket_name = c.split(".")
            if task_name == "context":
                conditions.append(self.ctx[socket_name])
            else:
                conditions.append(self.ctx.tasks[task_name]["results"][socket_name])
        print("conditions: ", conditions)
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_task_state(condition_tasks, "SKIPPED")
        return should_run

    def check_for_conditions(self) -> bool:
        print("Is a for workgraph")
        condition_tasks = [c[0] for c in self.ctx.workgraph["conditions"]]
        self.run_taskss(condition_tasks)
        conditions = [self.ctx._count < len(self.ctx.sequence)] + [
            self.ctx.tasks[c[0]]["results"][c[1]]
            for c in self.ctx.workgraph["conditions"]
        ]
        print("conditions: ", conditions)
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_task_state(condition_tasks, "SKIPPED")
            self.ctx["i"] = self.ctx.sequence[self.ctx._count]
        self.ctx._count += 1
        return should_run

    def run_taskss(self, names: t.List[str], continue_workgraph: bool = True) -> None:
        """Run task
        Here we use ToContext to pass the results of the run to the next step.
        This will force the engine to wait for all the submitted processes to
        finish before continuing to the next step.
        """
        from aiida_workgraph.utils import (
            get_executor,
            create_data_node,
            update_nested_dict_with_special_keys,
        )

        for name in names:
            print("-" * 60)
            task = self.ctx.tasks[name]
            if task["metadata"]["node_type"].upper() in [
                "CALCJOB",
                "WORKCHAIN",
                "GRAPH_BUILDER",
                "WORKGRAPH",
                "PYTHONTASK",
                "SHELLTASK",
            ]:
                if len(self._awaitables) > self.ctx.max_number_awaitables:
                    print(
                        MAX_NUMBER_AWAITABLES_MSG.format(
                            self.ctx.max_number_awaitables, name
                        )
                    )
                    continue
            self.report(f"Run task: {name}, type: {task['metadata']['node_type']}")
            # print("Run task: ", name)
            # print("executor: ", task["executor"])
            executor, _ = get_executor(task["executor"])
            print("executor: ", executor)
            args, kwargs, var_args, var_kwargs, args_dict = self.get_inputs(task)
            for i, key in enumerate(self.ctx.tasks[name]["metadata"]["args"]):
                kwargs[key] = args[i]
            # update the port namespace
            kwargs = update_nested_dict_with_special_keys(kwargs)
            print("args: ", args)
            print("kwargs: ", kwargs)
            print("var_kwargs: ", var_kwargs)
            # kwargs["meta.label"] = name
            # output must be a Data type or a mapping of {string: Data}
            task["results"] = {}
            if task["metadata"]["node_type"].upper() == "NODE":
                print("task  type: node.")
                results = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
                task["process"] = results
                task["results"] = {task["outputs"][0]["name"]: results}
                self.ctx.input_tasks[name] = results
                self.ctx.tasks[name]["state"] = "FINISHED"
                self.task_to_context(name)
                # ValueError: attempted to add an input link after the process node was already stored.
                # self.node.base.links.add_incoming(results, "INPUT_WORK", name)
                self.report(f"Task: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif task["metadata"]["node_type"].upper() == "DATA":
                print("task  type: data.")
                for key in self.ctx.tasks[name]["metadata"]["args"]:
                    kwargs.pop(key, None)
                results = create_data_node(executor, args, kwargs)
                task["results"] = {task["outputs"][0]["name"]: results}
                task["process"] = results
                self.ctx.new_data[name] = results
                self.ctx.tasks[name]["state"] = "FINISHED"
                self.task_to_context(name)
                self.report(f"Task: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif task["metadata"]["node_type"].upper() in [
                "CALCFUNCTION",
                "WORKFUNCTION",
            ]:
                print("task type: calcfunction/workfunction.")
                kwargs.setdefault("metadata", {})
                kwargs["metadata"].update({"call_link_label": name})
                try:
                    # since aiida 2.5.0, we need to use args_dict to pass the args to the run_get_node
                    if var_kwargs is None:
                        results, process = run_get_node(executor, **kwargs)
                    else:
                        results, process = run_get_node(
                            executor, **kwargs, **var_kwargs
                        )
                    process.label = name
                    # only one output
                    if isinstance(results, orm.Data):
                        results = {task["outputs"][0]["name"]: results}
                    task["results"] = results
                    # print("results: ", results)
                    task["process"] = process
                    self.ctx.tasks[name]["state"] = "FINISHED"
                    self.task_to_context(name)
                    self.report(f"Task: {name} finished.")
                except Exception as e:
                    print(e)
                    self.report(e)
                    self.ctx.tasks[name]["state"] = "FAILED"
                    # set child state to FAILED
                    self.set_task_state(
                        self.ctx.connectivity["child_node"][name], "SKIPPED"
                    )
                    print(f"Task: {name} failed.")
                    self.report(f"Task: {name} failed.")
                # exclude the current tasks from the next run
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif task["metadata"]["node_type"].upper() in ["CALCJOB", "WORKCHAIN"]:
                # process = run_get_node(executor, *args, **kwargs)
                print("task type: calcjob/workchain.")
                kwargs.setdefault("metadata", {})
                kwargs["metadata"].update({"call_link_label": name})
                # transfer the args to kwargs
                process = self.submit(executor, **kwargs)
                process.label = name
                task["process"] = process
                self.ctx.tasks[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif task["metadata"]["node_type"].upper() in ["GRAPH_BUILDER"]:
                print("task type: graph_builder.")
                wg = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
                wg.name = name
                wg.group_outputs = self.ctx.tasks[name]["metadata"]["group_outputs"]
                wg.parent_uuid = self.node.uuid
                wg.save(metadata={"call_link_label": name})
                print("submit workgraph: ")
                process = self.submit(wg.process_inited)
                task["process"] = process
                self.ctx.tasks[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif task["metadata"]["node_type"].upper() in ["WORKGRAPH"]:
                from .utils import prepare_for_workgraph_task
                from aiida_workgraph.utils.analysis import WorkGraphSaver

                inputs, wgdata = prepare_for_workgraph_task(task, kwargs)
                process_inited = WorkGraphEngine(inputs=inputs)
                process_inited.runner.persister.save_checkpoint(process_inited)
                saver = WorkGraphSaver(process_inited.node, wgdata)
                saver.save()
                print("submit workgraph: ")
                process = self.submit(process_inited)
                task["process"] = process
                self.ctx.tasks[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif task["metadata"]["node_type"].upper() in ["PYTHONTASK"]:
                from aiida_workgraph.calculations.python import PythonTask
                from .utils import prepare_for_python_task

                inputs = prepare_for_python_task(task, kwargs, var_kwargs)
                # since aiida 2.5.0, we can pass inputs directly to the submit, no need to use **inputs
                process = self.submit(
                    PythonTask,
                    **inputs,
                )
                process.label = name
                task["process"] = process
                self.ctx.tasks[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif task["metadata"]["node_type"].upper() in ["SHELLTASK"]:
                from aiida_shell.calculations.shell import ShellJob
                from .utils import prepare_for_shell_task

                inputs = prepare_for_shell_task(task, kwargs)
                # since aiida 2.5.0, we can pass inputs directly to the submit, no need to use **inputs
                process = self.submit(
                    ShellJob,
                    **inputs,
                )
                process.label = name
                task["process"] = process
                self.ctx.tasks[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif task["metadata"]["node_type"].upper() in ["NORMAL"]:
                print("Task  type: Normal.")
                # normal function does not have a process
                if "context" in task["metadata"]["kwargs"]:
                    self.ctx.task_name = name
                    kwargs.update({"context": self.ctx})
                for key in self.ctx.tasks[name]["metadata"]["args"]:
                    kwargs.pop(key, None)
                results = self.run_executor(
                    executor, args, kwargs, var_args, var_kwargs
                )
                # task["process"] = results
                if isinstance(results, tuple):
                    if len(task["outputs"]) != len(results):
                        return self.exit_codes.OUTPUS_NOT_MATCH_RESULTS
                    for i in range(len(task["outputs"])):
                        task["results"][task["outputs"][i]["name"]] = results[i]
                elif isinstance(results, dict):
                    task["results"] = results
                else:
                    task["results"][task["outputs"][0]["name"]] = results
                # save the results to the database (as a extra field of the task)
                # this is disabled
                # self.save_results_to_extras(name)
                self.ctx.input_tasks[name] = results
                self.ctx.tasks[name]["state"] = "FINISHED"
                self.task_to_context(name)
                self.report(f"Task: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
                # print("result from node: ", task["results"])
            else:
                print("Task type: unknown.")
                # self.report("Unknow task type {}".format(task["metadata"]["node_type"]))
                return self.exit_codes.UNKNOWN_TASK_TYPE

    def get_inputs(
        self, task: t.Dict[str, t.Any]
    ) -> t.Tuple[
        t.List[t.Any],
        t.Dict[str, t.Any],
        t.Optional[t.List[t.Any]],
        t.Optional[t.Dict[str, t.Any]],
        t.Dict[str, t.Any],
    ]:
        """Get input based on the links."""
        from aiida_workgraph.utils import get_nested_dict

        args = []
        args_dict = {}
        kwargs = {}
        var_args = None
        var_kwargs = None
        properties = task.get("properties", {})
        # TODO: check if input is linked, otherwise use the property value
        inputs = {}
        for input in task["inputs"]:
            # print(f"input: {input['name']}")
            if len(input["links"]) == 0:
                inputs[input["name"]] = self.update_context_variable(
                    properties[input["name"]]["value"]
                )
            elif len(input["links"]) == 1:
                link = input["links"][0]
                if self.ctx.tasks[link["from_node"]]["results"] is None:
                    inputs[input["name"]] = None
                else:
                    # handle the special socket _wait, _outputs
                    if link["from_socket"] == "_wait":
                        continue
                    elif link["from_socket"] == "_outputs":
                        inputs[input["name"]] = self.ctx.tasks[link["from_node"]][
                            "results"
                        ]
                    else:
                        inputs[input["name"]] = get_nested_dict(
                            self.ctx.tasks[link["from_node"]]["results"],
                            link["from_socket"],
                        )
            # handle the case of multiple outputs
            elif len(input["links"]) > 1:
                value = {}
                for link in input["links"]:
                    name = f'{link["from_node"]}_{link["from_socket"]}'
                    # handle the special socket _wait, _outputs
                    if link["from_socket"] == "_wait":
                        continue
                    if self.ctx.tasks[link["from_node"]]["results"] is None:
                        value[name] = None
                    else:
                        value[name] = self.ctx.tasks[link["from_node"]]["results"][
                            link["from_socket"]
                        ]
                inputs[input["name"]] = value
        for name in task["metadata"].get("args", []):
            if name in inputs:
                args.append(inputs[name])
                args_dict[name] = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                args.append(value)
                args_dict[name] = value
        for name in task["metadata"].get("kwargs", []):
            if name in inputs:
                kwargs[name] = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                kwargs[name] = value
        if task["metadata"]["var_args"] is not None:
            name = task["metadata"]["var_args"]
            if name in inputs:
                var_args = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                var_args = value
        if task["metadata"]["var_kwargs"] is not None:
            name = task["metadata"]["var_kwargs"]
            if name in inputs:
                var_kwargs = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                var_kwargs = value
        return args, kwargs, var_args, var_kwargs, args_dict

    def update_context_variable(self, value: t.Any) -> t.Any:
        # replace context variables
        from aiida_workgraph.utils import get_nested_dict

        """Get value from context."""
        if isinstance(value, dict):
            for key, sub_value in value.items():
                value[key] = self.update_context_variable(sub_value)
        elif (
            isinstance(value, str)
            and value.strip().startswith("{{")
            and value.strip().endswith("}}")
        ):
            name = value[2:-2].strip()
            return get_nested_dict(self.ctx, name)
        return value

    def task_to_context(self, name: str) -> None:
        """Export task result to context."""
        from aiida_workgraph.utils import update_nested_dict
        from aiida.common.exceptions import NotExistentKeyError

        items = self.ctx.tasks[name]["to_context"]
        for item in items:
            try:
                result = self.ctx.tasks[name]["results"][item[0]]
                update_nested_dict(self.ctx, item[1], result)
            except NotExistentKeyError as e:
                print(f"Warning: {e}. Skipping update for item {item[0]}")

    def check_task_state(self, name: str) -> None:
        """Check task states.

        - if all input tasks finished, launch task
        - if task is a scatter task, check if all scattered tasks finished
        """
        # print(f"    Check task {name} state: ")
        if self.ctx.tasks[name]["state"] in ["CREATED", "WAITING"]:
            ready, output = self.check_parent_state(name)
            if ready:
                # print(f"    Task {name} is ready to launch.")
                self.ctx.msgs.append(f"task,{name}:action:LAUNCH")  # noqa E231
        elif self.ctx.tasks[name]["state"] in ["SCATTERED"]:
            state, action = self.check_scattered_state(name)
            self.ctx.msgs.append(f"task,{name}:state:{state}")  # noqa E231
        else:
            # print(f"    Task {name} is in state {self.ctx.tasks[name]['state']}")
            pass

    def check_parent_state(self, name: str) -> t.Tuple[bool, t.Optional[str]]:
        task = self.ctx.tasks[name]
        inputs = task.get("inputs", None)
        wait_tasks = self.ctx.tasks[name].get("wait", [])
        # print("    wait_tasks: ", wait_tasks)
        ready = True
        if inputs is None and len(wait_tasks) == 0:
            return ready
        else:
            # check the wait task first
            for task_name in wait_tasks:
                # in case the task is removed
                if task_name not in self.ctx.tasks:
                    continue
                if self.ctx.tasks[task_name]["state"] not in [
                    "FINISHED",
                    "SKIPPED",
                    "FAILED",
                ]:
                    ready = False
                    return ready, f"Task {name} wait for {task_name}"
            for input in inputs:
                # print("    input, ", input["from_node"], self.ctx.tasks[input["from_node"]]["state"])
                for link in input["links"]:
                    if self.ctx.tasks[link["from_node"]]["state"] not in [
                        "FINISHED",
                        "SKIPPED",
                        "FAILED",
                    ]:
                        ready = False
                        return (
                            ready,
                            f"{name}, input: {link['from_node']} is {self.ctx.tasks[link['from_node']]['state']}",
                        )
        return ready, None

    # def expose_graph_build_outputs(self, name):
    #     # print("expose_graph_build_outputs")
    #     outputs = {}
    #     process = self.ctx.tasks[name]["process"]
    #     outgoing = process.base.links.get_outgoing()
    #     for output in self.ctx.tasks[name]["group_outputs"]:
    #         node = outgoing.get_node_by_label(output[0])
    #         outputs[output[2]] = getattr(node.outputs, output[1])
    #     return outputs
    def reset(self) -> None:
        print("Reset")
        self.ctx._execution_count += 1
        self.set_task_state(self.ctx.tasks.keys(), "CREATED")

    def set_task_state(
        self, names: t.Union[t.List[str], t.Sequence[str]], value: str
    ) -> None:
        """Set task state"""
        for name in names:
            self.ctx.tasks[name]["state"] = value

    def run_executor(
        self,
        executor: t.Callable,
        args: t.List[t.Any],
        kwargs: t.Dict[str, t.Any],
        var_args: t.Optional[t.List[t.Any]],
        var_kwargs: t.Optional[t.Dict[str, t.Any]],
    ) -> t.Any:
        if var_kwargs is None:
            return executor(*args, **kwargs)
        else:
            print("var_kwargs: ", var_kwargs)
            return executor(*args, **kwargs, **var_kwargs)

    def save_results_to_extras(self, name: str) -> None:
        """Save the results to the base.extras.
        For the outputs of a Normal task, they are not saved to the database like the calcjob or workchain.
        One temporary solution is to save the results to the base.extras. In order to do this, we need to
        serialize the results
        """
        from aiida_workgraph.utils import get_executor

        results = self.ctx.tasks[name]["results"]
        if results is None:
            return
        datas = {}
        for key, value in results.items():
            # find outptus sockets with the name as key
            output = [
                output
                for output in self.ctx.tasks[name]["outputs"]
                if output["name"] == key
            ]
            if len(output) == 0:
                continue
            output = output[0]
            Executor, _ = get_executor(output["serialize"])
            datas[key] = Executor(value)
        self.node.set_extra(f"nodes__results__{name}", datas)

    def finalize(self) -> t.Optional[ExitCode]:
        """"""
        from aiida_workgraph.utils import get_nested_dict, update_nested_dict

        # expose group outputs
        group_outputs = {}
        print("group outputs: ", self.ctx.workgraph["metadata"]["group_outputs"])
        for output in self.ctx.workgraph["metadata"]["group_outputs"]:
            print("output: ", output)
            task_name, socket_name = output[0].split(".")
            if task_name == "context":
                update_nested_dict(
                    group_outputs, output[1], get_nested_dict(self.ctx, socket_name)
                )
            else:
                update_nested_dict(
                    group_outputs,
                    output[1],
                    self.ctx.tasks[task_name]["results"][socket_name],
                )
        self.out("group_outputs", group_outputs)
        self.out("new_data", self.ctx.new_data)
        self.out("execution_count", orm.Int(self.ctx._execution_count).store())
        self.report("Finalize")
        for name, task in self.ctx.tasks.items():
            if task["state"] == "FAILED":
                print(f"    Task {name} failed.")
                return self.exit_codes.TASK_FAILED
        print(f"Finalize workgraph {self.ctx.workgraph['name']}\n")
