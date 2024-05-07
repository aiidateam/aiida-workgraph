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
from aiida.orm import Node, ProcessNode, WorkChainNode
from aiida.orm.utils import load_node


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
class WorkGraph(Process, metaclass=Protect):
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
        spec.input_namespace("input_nodes", dynamic=True, required=False)
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
        spec.exit_code(202, "UNKNOWN_NODE_TYPE", message="The node type is unknown.")
        #
        spec.exit_code(
            301,
            "OUTPUS_NOT_MATCH_RESULTS",
            message="The outputs of the process do not match the results.",
        )
        spec.exit_code(
            302,
            "NODE_FAILED",
            message="Some of the nodes failed.",
        )
        spec.exit_code(
            303,
            "NODE_NON_ZERO_EXIT_STATUS",
            message="Some of the nodes exited with non-zero status.",
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
            finished = self.is_workgraph_finished()

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

        # node finished, update the node state and result
        # udpate the node state
        self.update_node_state(awaitable.key)
        # try to resume the workgraph, if the workgraph is already resumed
        # by other awaitable, this will not work
        try:
            self.resume()
        except Exception as e:
            print(e)

    def setup(self) -> None:
        # track if the awaitable callback is added to the runner
        self.ctx._awaitable_actions = []
        self.ctx.new_data = dict()
        self.ctx.input_nodes = dict()
        # read the latest workgraph data
        wgdata = self.read_wgdata_from_base()
        self.init_ctx(wgdata)
        #
        self.ctx.msgs = []
        self.node.set_process_label(f"WorkGraph: {self.ctx.workgraph['name']}")
        self.ctx._execution_count = 0
        # init node results
        self.set_node_results()
        # while workgraph
        if self.ctx.workgraph["workgraph_type"].upper() == "WHILE":
            self.ctx._max_iteration = self.ctx.workgraph.get("max_iteration", 1000)
            should_run = self.check_while_conditions()
            if not should_run:
                self.set_node_state(self.ctx.nodes.keys(), "SKIPPED")
        # for workgraph
        if self.ctx.workgraph["workgraph_type"].upper() == "FOR":
            should_run = self.check_for_conditions()
            if not should_run:
                self.set_node_state(self.ctx.nodes.keys(), "SKIPPED")

    def setup_ctx_workgraph(self, wgdata: t.Dict[str, t.Any]) -> None:
        """setup the workgraph in the context."""
        self.ctx.nodes = wgdata["nodes"]
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
        for name, node in wgdata["nodes"].items():
            node["state"] = self.ctx.nodes[name]["state"]
            node["results"] = self.ctx.nodes[name].get("results")
            node["process"] = self.ctx.nodes[name].get("process")
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

    def set_node_results(self) -> None:
        for _, node in self.ctx.nodes.items():
            if node.get("process"):
                if isinstance(node["process"], str):
                    node["process"] = orm.load_node(node["process"])
                self.set_node_result(node)
            self.set_node_result(node)

    def set_node_result(self, node: t.Dict[str, t.Any]) -> None:
        name = node["name"]
        # print(f"set node result: {name}")
        if node.get("process"):
            # print(f"set node result: {name} process")
            state = node["process"].process_state.value.upper()
            if state == "FINISHED":
                node["state"] = state
                if node["metadata"]["node_type"] == "graph_builder":
                    # expose the outputs of workgraph
                    node["results"] = getattr(
                        node["process"].outputs, "group_outputs", None
                    )
                    # self.ctx.new_data[name] = outputs
                elif node["metadata"]["node_type"] == "workgraph":
                    # expose the outputs of all the nodes in the workgraph
                    node["results"] = {}
                    outgoing = node["process"].base.links.get_outgoing()
                    for link in outgoing.all():
                        if isinstance(link.node, ProcessNode) and getattr(
                            link.node, "process_state", False
                        ):
                            node["results"][link.link_label] = link.node.outputs
                else:
                    node["results"] = node["process"].outputs
                    # self.ctx.new_data[name] = node["results"]
                self.ctx.nodes[name]["state"] = "FINISHED"
                self.node_to_context(name)
                self.report(f"Node: {name} finished.")
            elif state == "EXCEPTED":
                node["state"] = state
                node["results"] = node["process"].outputs
                # self.ctx.new_data[name] = node["results"]
                self.ctx.nodes[name]["state"] = "FAILED"
                # set child state to FAILED
                self.set_node_state(self.ctx.connectivity["child_node"][name], "FAILED")
                self.report(f"Node: {name} failed.")
        else:
            node["results"] = None

    def apply_actions(self) -> None:
        """Apply actions to the workgraph.
        The actions are stored in the base.extras["workgraph_queue"].
        The index of the last applied action is stored in the base.extras["workgraph_queue_index"].
        """
        msgs = self.node.base.extras.get("workgraph_queue", [])
        index = self.node.base.extras.get("workgraph_queue_index", 0)
        for msg in msgs[index:]:
            header, msg = msg.split(",")
            if header == "node":
                self.apply_node_actions(msg)
            else:
                self.report(f"Unknow message type {msg}")
            index += 1
            self.report("Apply actions: {}".format(msg))
            msgs = self.node.base.extras.set("workgraph_queue_index", index)

    def apply_node_actions(self, msg: str) -> None:
        """Apply node actions to the workgraph."""
        name, action = msg.split(":")
        print("apply node actions: ", name, action)
        if action == "RESET":
            self.reset_node(name)

    def reset_node(self, name: str) -> None:
        """Reset node."""
        self.ctx.nodes[name]["state"] = "CREATED"
        # reset its child nodes
        names = self.ctx.connectivity["child_node"][name]
        for name in names:
            self.ctx.nodes[name]["state"] = "CREATED"
            self.ctx.nodes[name]["result"] = None
            self.ctx.nodes[name]["process"] = None

    def continue_workgraph(self, exclude: t.Optional[t.List[str]] = None) -> None:
        print("Continue workgraph.")
        exclude = exclude or []
        self.report("Continue workgraph.")
        self.update_workgraph_from_base()
        self.apply_actions()
        node_to_run = []
        for name, node in self.ctx.nodes.items():
            # update node state
            if node["state"] in ["RUNNING", "FINISHED", "FAILED", "SKIPPED"]:
                continue
            ready, output = self.check_parent_state(name)
            if ready and name not in exclude:
                node_to_run.append(name)
        #
        self.report("nodes ready to run: {}".format(",".join(node_to_run)))
        self.run_nodes(node_to_run)

    def update_node_state(self, name: str) -> None:
        """Update ndoe state if node is a Awaitable."""
        print("update node state: ", name)
        node = self.ctx.nodes[name]
        if (
            node["metadata"]["node_type"]
            in [
                "calcfunction",
                "workfunction",
                "calcjob",
                "workchain",
                "graph_builder",
                "workgraph",
            ]
            and node["state"] == "RUNNING"
        ):
            self.set_node_result(node)

    def is_workgraph_finished(self) -> bool:
        """Check if the workgraph is finished.
        For `while` workgraph, we need check its conditions"""
        is_finished = True
        for name, node in self.ctx.nodes.items():
            # self.update_node_state(name)
            print("node: ", name, node["state"])
            if node["state"] in ["RUNNING", "CREATED", "READY"]:
                is_finished = False
        if is_finished:
            if self.ctx.workgraph["workgraph_type"].upper() == "WHILE":
                should_run = self.check_while_conditions()
                is_finished = not should_run
            if self.ctx.workgraph["workgraph_type"].upper() == "FOR":
                should_run = self.check_for_conditions()
                is_finished = not should_run
        print("is workgraph finished: ", is_finished)
        return is_finished

    def check_while_conditions(self) -> bool:
        """Check while conditions.
        Run all condition nodes and check if all the conditions are True.
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
        condition_nodes = []
        for c in self.ctx.workgraph["conditions"]:
            node_name, socket_name = c.split(".")
            if "node_name" != "context":
                condition_nodes.append(node_name)
        self.run_nodes(condition_nodes, continue_workgraph=False)
        conditions = []
        for c in self.ctx.workgraph["conditions"]:
            node_name, socket_name = c.split(".")
            if node_name == "context":
                conditions.append(self.ctx[socket_name])
            else:
                conditions.append(self.ctx.nodes[node_name]["results"][socket_name])
        print("conditions: ", conditions)
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_node_state(condition_nodes, "SKIPPED")
        return should_run

    def check_for_conditions(self) -> bool:
        print("Is a for workgraph")
        condition_nodes = [c[0] for c in self.ctx.workgraph["conditions"]]
        self.run_nodes(condition_nodes)
        conditions = [self.ctx._count < len(self.ctx.sequence)] + [
            self.ctx.nodes[c[0]]["results"][c[1]]
            for c in self.ctx.workgraph["conditions"]
        ]
        print("conditions: ", conditions)
        should_run = False not in conditions
        if should_run:
            self.reset()
            self.set_node_state(condition_nodes, "SKIPPED")
            self.ctx["i"] = self.ctx.sequence[self.ctx._count]
        self.ctx._count += 1
        return should_run

    def run_nodes(self, names: t.List[str], continue_workgraph: bool = True) -> None:
        """Run node
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
            node = self.ctx.nodes[name]
            if node["metadata"]["node_type"] in [
                "calcjob",
                "workchain",
                "graph_builder",
                "workgraph",
            ]:
                if len(self._awaitables) > self.ctx.max_number_awaitables:
                    print(
                        MAX_NUMBER_AWAITABLES_MSG.format(
                            self.ctx.max_number_awaitables, name
                        )
                    )
                    continue
            self.report(f"Run node: {name}, type: {node['metadata']['node_type']}")
            # print("Run node: ", name)
            # print("executor: ", node["executor"])
            executor, _ = get_executor(node["executor"])
            print("executor: ", executor)
            args, kwargs, var_args, var_kwargs, args_dict = self.get_inputs(node)
            for i, key in enumerate(self.ctx.nodes[name]["metadata"]["args"]):
                kwargs[key] = args[i]
            # update the port namespace
            kwargs = update_nested_dict_with_special_keys(kwargs)
            # print("args: ", args)
            # print("kwargs: ", kwargs)
            # print("var_kwargs: ", var_kwargs)
            # kwargs["meta.label"] = name
            # output must be a Data type or a mapping of {string: Data}
            node["results"] = {}
            if node["metadata"]["node_type"] == "node":
                print("node  type: node.")
                results = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
                node["process"] = results
                node["results"] = {node["outputs"][0]["name"]: results}
                self.ctx.input_nodes[name] = results
                self.ctx.nodes[name]["state"] = "FINISHED"
                self.node_to_context(name)
                # ValueError: attempted to add an input link after the process node was already stored.
                # self.node.base.links.add_incoming(results, "INPUT_WORK", name)
                self.report(f"Node: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif node["metadata"]["node_type"] == "data":
                print("node  type: data.")
                for key in self.ctx.nodes[name]["metadata"]["args"]:
                    kwargs.pop(key, None)
                results = create_data_node(executor, args, kwargs)
                node["results"] = {node["outputs"][0]["name"]: results}
                node["process"] = results
                self.ctx.new_data[name] = results
                self.ctx.nodes[name]["state"] = "FINISHED"
                self.node_to_context(name)
                self.report(f"Node: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif node["metadata"]["node_type"] in ["calcfunction", "workfunction"]:
                print("node type: calcfunction/workfunction.")
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
                        results = {node["outputs"][0]["name"]: results}
                    node["results"] = results
                    # print("results: ", results)
                    node["process"] = process
                    self.ctx.nodes[name]["state"] = "FINISHED"
                    self.node_to_context(name)
                    self.report(f"Node: {name} finished.")
                except Exception as e:
                    print(e)
                    self.report(e)
                    self.ctx.nodes[name]["state"] = "FAILED"
                    # set child state to FAILED
                    self.set_node_state(
                        self.ctx.connectivity["child_node"][name], "FAILED"
                    )
                    print(f"Node: {name} failed.")
                    self.report(f"Node: {name} failed.")
                # exclude the current nodes from the next run
                if continue_workgraph:
                    self.continue_workgraph(names)
            elif node["metadata"]["node_type"] in ["calcjob", "workchain"]:
                # process = run_get_node(executor, *args, **kwargs)
                print("node  type: calcjob/workchain.")
                kwargs.setdefault("metadata", {})
                kwargs["metadata"].update({"call_link_label": name})
                # transfer the args to kwargs
                process = self.submit(executor, **kwargs)
                process.label = name
                node["process"] = process
                self.ctx.nodes[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif node["metadata"]["node_type"] in ["graph_builder"]:
                print("node  type: graph_builder.")
                wg = self.run_executor(executor, [], kwargs, var_args, var_kwargs)
                wg.name = name
                wg.group_outputs = self.ctx.nodes[name]["metadata"]["group_outputs"]
                wg.parent_uuid = self.node.uuid
                wg.save(metadata={"call_link_label": name})
                print("submit workgraph: ")
                process = self.submit(wg.process_inited)
                node["process"] = process
                self.ctx.nodes[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif node["metadata"]["node_type"] in ["workgraph"]:
                from aiida_workgraph.utils import merge_properties
                from aiida_workgraph.utils.analysis import WorkGraphSaver

                print("node  type: workgraph.")
                wgdata = node["executor"]["wgdata"]
                wgdata["name"] = name
                wgdata["metadata"]["group_outputs"] = self.ctx.nodes[name]["metadata"][
                    "group_outputs"
                ]
                # update the workgraph data by kwargs
                for node_name, data in kwargs.items():
                    # because kwargs is updated using update_nested_dict_with_special_keys
                    # which means the data is grouped by the node name
                    for socket_name, value in data.items():
                        wgdata["nodes"][node_name]["properties"][socket_name][
                            "value"
                        ] = value
                # merge the properties
                merge_properties(wgdata)
                metadata = {"call_link_label": name}
                inputs = {"wg": wgdata, "metadata": metadata}
                process_inited = WorkGraph(inputs=inputs)
                process_inited.runner.persister.save_checkpoint(process_inited)
                saver = WorkGraphSaver(process_inited.node, wgdata)
                saver.save()
                print("submit workgraph: ")
                process = self.submit(process_inited)
                node["process"] = process
                self.ctx.nodes[name]["state"] = "RUNNING"
                self.to_context(**{name: process})
            elif node["metadata"]["node_type"] in ["Normal"]:
                print("node  type: Normal.")
                # normal function does not have a process
                if "context" in node["metadata"]["kwargs"]:
                    self.ctx.node_name = name
                    kwargs.update({"context": self.ctx})
                for key in self.ctx.nodes[name]["metadata"]["args"]:
                    kwargs.pop(key, None)
                results = self.run_executor(
                    executor, args, kwargs, var_args, var_kwargs
                )
                # node["process"] = results
                if isinstance(results, tuple):
                    if len(node["outputs"]) != len(results):
                        return self.exit_codes.OUTPUS_NOT_MATCH_RESULTS
                    for i in range(len(node["outputs"])):
                        node["results"][node["outputs"][i]["name"]] = results[i]
                elif isinstance(results, dict):
                    node["results"] = results
                else:
                    node["results"][node["outputs"][0]["name"]] = results
                # save the results to the database (as a extra field of the node)
                # this is disabled
                # self.save_results_to_extras(name)
                self.ctx.input_nodes[name] = results
                self.ctx.nodes[name]["state"] = "FINISHED"
                self.node_to_context(name)
                self.report(f"Node: {name} finished.")
                if continue_workgraph:
                    self.continue_workgraph(names)
                # print("result from node: ", node["results"])
            else:
                print("node  type: unknown.")
                # self.report("Unknow node type {}".format(node["metadata"]["node_type"]))
                return self.exit_codes.UNKNOWN_NODE_TYPE

    def get_inputs(
        self, node: t.Dict[str, t.Any]
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
        properties = node.get("properties", {})
        # TODO: check if input is linked, otherwise use the property value
        inputs = {}
        for input in node["inputs"]:
            # print(f"input: {input['name']}")
            if len(input["links"]) == 0:
                inputs[input["name"]] = self.update_context_variable(
                    properties[input["name"]]["value"]
                )
            elif len(input["links"]) == 1:
                link = input["links"][0]
                if self.ctx.nodes[link["from_node"]]["results"] is None:
                    inputs[input["name"]] = None
                else:
                    # handle the special socket _wait, _outputs
                    if link["from_socket"] == "_wait":
                        continue
                    elif link["from_socket"] == "_outputs":
                        inputs[input["name"]] = self.ctx.nodes[link["from_node"]][
                            "results"
                        ]
                    else:
                        inputs[input["name"]] = get_nested_dict(
                            self.ctx.nodes[link["from_node"]]["results"],
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
                    if self.ctx.nodes[link["from_node"]]["results"] is None:
                        value[name] = None
                    else:
                        value[name] = self.ctx.nodes[link["from_node"]]["results"][
                            link["from_socket"]
                        ]
                inputs[input["name"]] = value
        for name in node["metadata"].get("args", []):
            if name in inputs:
                args.append(inputs[name])
                args_dict[name] = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                args.append(value)
                args_dict[name] = value
        for name in node["metadata"].get("kwargs", []):
            if name in inputs:
                kwargs[name] = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                kwargs[name] = value
        if node["metadata"]["var_args"] is not None:
            name = node["metadata"]["var_args"]
            if name in inputs:
                var_args = inputs[name]
            else:
                value = self.update_context_variable(properties[name]["value"])
                var_args = value
        if node["metadata"]["var_kwargs"] is not None:
            name = node["metadata"]["var_kwargs"]
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

    def node_to_context(self, name: str) -> None:
        """Export node result to context."""
        from aiida_workgraph.utils import update_nested_dict

        items = self.ctx.nodes[name]["to_context"]
        for item in items:
            update_nested_dict(
                self.ctx, item[1], self.ctx.nodes[name]["results"][item[0]]
            )

    def check_node_state(self, name: str) -> None:
        """Check node states.

        - if all input nodes finished, launch node
        - if node is a scatter node, check if all scattered nodes finished
        """
        # print(f"    Check node {name} state: ")
        if self.ctx.nodes[name]["state"] in ["CREATED", "WAITING"]:
            ready, output = self.check_parent_state(name)
            if ready:
                # print(f"    Node {name} is ready to launch.")
                self.ctx.msgs.append(f"node,{name}:action:LAUNCH")  # noqa E231
        elif self.ctx.nodes[name]["state"] in ["SCATTERED"]:
            state, action = self.check_scattered_state(name)
            self.ctx.msgs.append(f"node,{name}:state:{state}")  # noqa E231
        else:
            # print(f"    Node {name} is in state {self.ctx.nodes[name]['state']}")
            pass

    def check_parent_state(self, name: str) -> t.Tuple[bool, t.Optional[str]]:
        node = self.ctx.nodes[name]
        inputs = node.get("inputs", None)
        wait_nodes = self.ctx.nodes[name].get("wait", [])
        # print("    wait_nodes: ", wait_nodes)
        ready = True
        if inputs is None and len(wait_nodes) == 0:
            return ready
        else:
            # check the wait node first
            for node_name in wait_nodes:
                # in case the node is removed
                if node_name not in self.ctx.nodes:
                    continue
                if self.ctx.nodes[node_name]["state"] not in [
                    "FINISHED",
                    "SKIPPED",
                    "FAILED",
                ]:
                    ready = False
                    return ready, f"Node {name} wait for {node_name}"
            for input in inputs:
                # print("    input, ", input["from_node"], self.ctx.nodes[input["from_node"]]["state"])
                for link in input["links"]:
                    if self.ctx.nodes[link["from_node"]]["state"] not in [
                        "FINISHED",
                        "SKIPPED",
                        "FAILED",
                    ]:
                        ready = False
                        return (
                            ready,
                            f"{name}, input: {link['from_node']} is {self.ctx.nodes[link['from_node']]['state']}",
                        )
        return ready, None

    # def expose_graph_build_outputs(self, name):
    #     # print("expose_graph_build_outputs")
    #     outputs = {}
    #     process = self.ctx.nodes[name]["process"]
    #     outgoing = process.base.links.get_outgoing()
    #     for output in self.ctx.nodes[name]["group_outputs"]:
    #         node = outgoing.get_node_by_label(output[0])
    #         outputs[output[2]] = getattr(node.outputs, output[1])
    #     return outputs
    def reset(self) -> None:
        print("Reset")
        self.ctx._execution_count += 1
        self.set_node_state(self.ctx.nodes.keys(), "CREATED")

    def set_node_state(
        self, names: t.Union[t.List[str], t.Sequence[str]], value: str
    ) -> None:
        """Set node state"""
        for name in names:
            self.ctx.nodes[name]["state"] = value

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
        For the outputs of a Normal node, they are not saved to the database like the calcjob or workchain.
        One temporary solution is to save the results to the base.extras. In order to do this, we need to
        serialize the results
        """
        from aiida_workgraph.utils import get_executor

        results = self.ctx.nodes[name]["results"]
        if results is None:
            return
        datas = {}
        for key, value in results.items():
            # find outptus sockets with the name as key
            output = [
                output
                for output in self.ctx.nodes[name]["outputs"]
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
            node_name, socket_name = output[0].split(".")
            if node_name == "context":
                update_nested_dict(
                    group_outputs, output[1], get_nested_dict(self.ctx, socket_name)
                )
            else:
                update_nested_dict(
                    group_outputs,
                    output[1],
                    self.ctx.nodes[node_name]["results"][socket_name],
                )
        self.out("group_outputs", group_outputs)
        self.out("new_data", self.ctx.new_data)
        self.out("execution_count", orm.Int(self.ctx._execution_count).store())
        self.report("Finalize")
        for name, node in self.ctx.nodes.items():
            if node["state"] == "FAILED":
                print(f"    Node {name} failed.")
                return self.exit_codes.NODE_FAILED
        print(f"Finalize workgraph {self.ctx.workgraph['name']}\n")
        # check if all nodes are finished with nonzero exit code
