"""AiiDA workflow components: WorkGraph."""
from __future__ import annotations

import collections.abc
import logging
import typing as t

from plumpy import process_comms
from plumpy.persistence import auto_persist
from plumpy.process_states import Continue, Wait
from plumpy.workchains import _PropagateReturn
import kiwipy

from aiida.common.extendeddicts import AttributeDict
from aiida.common.lang import override
from aiida import orm
from aiida.orm import Node, WorkChainNode
from aiida.orm.utils.serialize import deserialize_unsafe

from aiida.engine.processes.exit_code import ExitCode
from aiida.engine.processes.process import Process
from aiida.engine.processes.workchains.workchain import Protect, WorkChainSpec
from aiida_workgraph.utils import get_nested_dict, update_nested_dict
from .context_manager import ContextManager
from .awaitable_manager import AwaitableManager
from .task_manager import TaskManager
from .error_handler_manager import ErrorHandlerManager
from aiida.engine.processes.workchains.awaitable import Awaitable

if t.TYPE_CHECKING:
    from aiida.engine.runners import Runner  # pylint: disable=unused-import

__all__ = "WorkGraph"


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
        self.ctx_manager = ContextManager(
            self._context, process=self, logger=self.logger
        )
        self.awaitable_manager = AwaitableManager(
            self._awaitables, self.runner, self.logger, self, self.ctx_manager
        )
        self.task_manager = TaskManager(
            self.ctx_manager, self.logger, self.runner, self, self.awaitable_manager
        )
        self.error_handler_manager = ErrorHandlerManager(
            self, self.ctx_manager, self.logger
        )

    @classmethod
    def define(cls, spec: WorkChainSpec) -> None:
        super().define(spec)
        spec.input("input_file", valid_type=orm.SinglefileData, required=False)
        spec.input_namespace(
            "wg", dynamic=True, required=False, help="WorkGraph inputs"
        )
        spec.input_namespace("input_tasks", dynamic=True, required=False)
        spec.exit_code(2, "ERROR_SUBPROCESS", message="A subprocess has failed.")

        spec.outputs.dynamic = True

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
        # TODO I don't know why we need to reinitialize the context, awaitables, and task_manager
        # Need to initialize the context, awaitables, and task_manager
        self.ctx_manager = ContextManager(
            self._context, process=self, logger=self.logger
        )
        self.awaitable_manager = AwaitableManager(
            self._awaitables, self.runner, self.logger, self, self.ctx_manager
        )
        self.task_manager = TaskManager(
            self.ctx_manager, self.logger, self.runner, self, self.awaitable_manager
        )
        self.error_handler_manager = ErrorHandlerManager(
            self, self.ctx_manager, self.logger
        )
        # "_awaitables" is auto persisted.
        if self._awaitables:
            # For the "ascyncio.tasks.Task" awaitable, because there are only in-memory,
            # we need to reset the tasks and so that they can be re-run again.
            should_resume = False
            for awaitable in self._awaitables:
                if awaitable.target == "asyncio.tasks.Task":
                    self.awaitable_manager.resolve_awaitable(awaitable, None)
                    self.report(f"reset awaitable task: {awaitable.key}")
                    self.task_manager.reset_task(awaitable.key)
                    should_resume = True
            if should_resume:
                self.awaitable_manager.update_process_status()
                self.resume()
            # For other awaitables, because they exist in the db, we only need to re-register the callbacks
            self.ctx._awaitable_actions = []
            self.awaitable_manager.action_awaitables()

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
            self.task_manager.continue_workgraph()
        except _PropagateReturn as exception:
            finished, result = True, exception.exit_code
        else:
            finished, result = self.task_manager.is_workgraph_finished()

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
            self.awaitable_manager.action_awaitables()
            self.report("Process status: {}".format(self.node.process_status))
        else:
            self.call_soon(self.resume)

    def _build_process_label(self) -> str:
        """Use the workgraph name as the process label."""
        return f"WorkGraph<{self.inputs.wg['name']}>"

    def on_create(self) -> None:
        """Called when a Process is created."""
        from aiida_workgraph.utils.analysis import WorkGraphSaver

        super().on_create()
        wgdata = self.inputs.wg._dict
        restart_process = (
            orm.load_node(wgdata["restart_process"].value)
            if wgdata.get("restart_process")
            else None
        )
        saver = WorkGraphSaver(self.node, wgdata, restart_process=restart_process)
        saver.save()
        self.node.label = wgdata["name"]

    def setup(self) -> None:
        """Setup the variables in the context."""
        # track if the awaitable callback is added to the runner
        self.ctx._awaitable_actions = []
        self.ctx._new_data = {}
        self.ctx._executed_tasks = []
        # read the latest workgraph data
        wgdata = self.read_wgdata_from_base()
        self.init_ctx(wgdata)
        #
        self.ctx._msgs = []
        self.ctx._execution_count = 1
        # init task results
        self.task_manager.set_task_results()
        # while workgraph
        if self.ctx._workgraph["workgraph_type"].upper() == "WHILE":
            self.ctx._max_iteration = self.ctx._workgraph.get("max_iteration", 1000)
            should_run = self.task_manager.check_while_conditions()
            if not should_run:
                self.task_manager.set_tasks_state(self.ctx._tasks.keys(), "SKIPPED")
        # for workgraph
        if self.ctx._workgraph["workgraph_type"].upper() == "FOR":
            should_run = self.task_manager.check_for_conditions()
            if not should_run:
                self.task_manager.set_tasks_state(self.ctx._tasks.keys(), "SKIPPED")

    def setup_ctx_workgraph(self, wgdata: t.Dict[str, t.Any]) -> None:
        """setup the workgraph in the context."""

        self.ctx._tasks = wgdata["tasks"]
        self.ctx._links = wgdata["links"]
        self.ctx._connectivity = wgdata["connectivity"]
        self.ctx._ctrl_links = wgdata["ctrl_links"]
        self.ctx._workgraph = wgdata
        self.ctx._error_handlers = wgdata["error_handlers"]

    def read_wgdata_from_base(self) -> t.Dict[str, t.Any]:
        """Read workgraph data from base.extras."""
        from aiida_workgraph.orm.function_data import PickledLocalFunction

        wgdata = self.node.base.extras.get("_workgraph")
        for name, task in wgdata["tasks"].items():
            wgdata["tasks"][name] = deserialize_unsafe(task)
            for _, input in wgdata["tasks"][name]["inputs"].items():
                if input["property"] is None:
                    continue
                prop = input["property"]
                if isinstance(prop["value"], PickledLocalFunction):
                    prop["value"] = prop["value"].value
        wgdata["error_handlers"] = deserialize_unsafe(wgdata["error_handlers"])
        wgdata["context"] = deserialize_unsafe(wgdata["context"])
        return wgdata

    def init_ctx(self, wgdata: t.Dict[str, t.Any]) -> None:
        """Init the context from the workgraph data."""
        from aiida_workgraph.utils import update_nested_dict

        # set up the context variables
        self.ctx._max_number_awaitables = (
            wgdata["max_number_jobs"] if wgdata["max_number_jobs"] else 1000000
        )
        self.ctx["_count"] = 0
        for key, value in wgdata["context"].items():
            key = key.replace("__", ".")
            update_nested_dict(self.ctx, key, value)
        # set up the workgraph
        self.setup_ctx_workgraph(wgdata)

    def apply_action(self, msg: dict) -> None:

        if msg["catalog"] == "task":
            self.task_manager.apply_task_actions(msg)
        else:
            self.report(f"Unknow message type {msg}")

    def message_receive(
        self, _comm: kiwipy.Communicator, msg: t.Dict[str, t.Any]
    ) -> t.Any:
        """
        Coroutine called when the process receives a message from the communicator

        :param _comm: the communicator that sent the message
        :param msg: the message
        :return: the outcome of processing the message, the return value will be sent back as a response to the sender
        """
        self.logger.debug(
            "Process<%s>: received RPC message with communicator '%s': %r",
            self.pid,
            _comm,
            msg,
        )

        intent = msg[process_comms.INTENT_KEY]

        if intent == process_comms.Intent.PLAY:
            return self._schedule_rpc(self.play)
        if intent == process_comms.Intent.PAUSE:
            return self._schedule_rpc(
                self.pause, msg=msg.get(process_comms.MESSAGE_KEY, None)
            )
        if intent == process_comms.Intent.KILL:
            return self._schedule_rpc(
                self.kill, msg=msg.get(process_comms.MESSAGE_KEY, None)
            )
        if intent == process_comms.Intent.STATUS:
            status_info: t.Dict[str, t.Any] = {}
            self.get_status_info(status_info)
            return status_info
        if intent == "custom":
            return self._schedule_rpc(self.apply_action, msg=msg)

        # Didn't match any known intents
        raise RuntimeError("Unknown intent")

    def finalize(self) -> t.Optional[ExitCode]:
        """"""
        # expose outputs of the workgraph
        group_outputs = {}
        for output in self.ctx._workgraph["metadata"]["group_outputs"]:
            names = output["from"].split(".", 1)
            if names[0] == "context":
                if len(names) == 1:
                    raise ValueError("The output name should be context.key")
                update_nested_dict(
                    group_outputs,
                    output["name"],
                    get_nested_dict(self.ctx, names[1]),
                )
            else:
                # expose the whole outputs of the tasks
                if len(names) == 1:
                    update_nested_dict(
                        group_outputs,
                        output["name"],
                        self.ctx._tasks[names[0]]["results"],
                    )
                else:
                    # expose one output of the task
                    # note, the output may not exist
                    if names[1] in self.ctx._tasks[names[0]]["results"]:
                        update_nested_dict(
                            group_outputs,
                            output["name"],
                            self.ctx._tasks[names[0]]["results"][names[1]],
                        )
        self.out_many(group_outputs)
        # output the new data
        if self.ctx._new_data:
            self.out("new_data", self.ctx._new_data)
        self.out("execution_count", orm.Int(self.ctx._execution_count).store())
        self.report("Finalize workgraph.")
        for _, task in self.ctx._tasks.items():
            if self.task_manager.get_task_state_info(task["name"], "state") == "FAILED":
                return self.exit_codes.TASK_FAILED
