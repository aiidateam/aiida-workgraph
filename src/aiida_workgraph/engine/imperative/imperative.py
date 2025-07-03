"""AiiDA workflow components: WorkGraph."""
from __future__ import annotations

import typing as t

from plumpy.persistence import auto_persist
from plumpy.process_states import Continue, Wait
from plumpy.workchains import _PropagateReturn
from aiida.engine.processes.exit_code import ExitCode
from aiida_workgraph.engine.workgraph import WorkGraphEngine, WorkGraphSpec
from aiida_workgraph.socket import TaskSocketNamespace, TaskSocket

__all__ = "WorkGraph"


@auto_persist("_awaitables")
class WorkGraphImperativeEngine(WorkGraphEngine):
    """The `WorkGraph` class is used to construct workflows in AiiDA."""

    def on_create(self) -> None:
        """Called when a Process is created."""
        from aiida_workgraph.utils.analysis import WorkGraphSaver
        from aiida_workgraph import WorkGraph

        super(WorkGraphEngine, self).on_create()
        name = self._raw_inputs[WorkGraphSpec.WORKGRAPH_DATA_KEY]["name"]
        flow = self.inputs[WorkGraphSpec.WORKGRAPH_DATA_KEY]["flow"]
        self.node.label = name
        wg = WorkGraph(name)
        wg.add_task(flow, name="_flow")
        saver = WorkGraphSaver(
            self.node,
            wg.prepare_inputs()[WorkGraphSpec.WORKGRAPH_DATA_KEY],
        )
        saver.save()

    def execute_flow(self) -> t.Any:
        from node_graph.executor import NodeExecutor
        import asyncio

        name = "_flow"
        task = self.wg.tasks[name]
        executor = NodeExecutor(**task.get_executor()).executor
        inputs = self.inputs.workgraph_data["function_inputs"]

        awaitable_target = asyncio.ensure_future(
            executor(**inputs),
            loop=self.loop,
        )
        awaitable = self.awaitable_manager.construct_awaitable_function(
            name, awaitable_target
        )
        self.task_manager.state_manager.set_task_runtime_info(name, "state", "RUNNING")
        # save the awaitable to the temp, so that we can kill it if needed
        self.awaitable_manager.not_persisted_awaitables[name] = awaitable_target
        self.awaitable_manager.to_context(**{name: awaitable})

    def setup(self) -> None:
        """Setup the workgraph engine."""
        super().setup()

        # Execute the flow
        self.execute_flow()

    def _do_step(self) -> t.Any:
        """Execute the next step in the workgraph and return the result.

        If any awaitables were created, the process will enter in the Wait state,
        otherwise it will go to Continue.
        """
        result: t.Any = None

        try:
            self.task_manager.continue_workgraph()
        except _PropagateReturn as exception:
            finished, result = True, exception.exit_code
        else:
            finished, result = self.task_manager.is_workgraph_finished()

        # If the workgraph is finished or the result is an ExitCode, we exit by returning
        if finished and len(self._awaitables) == 0:
            if isinstance(result, ExitCode):
                return result
            else:
                return self.finalize()

        if self._awaitables:
            self.awaitable_manager.action_awaitables()
            return Wait(self._do_step, "Waiting before next step")

        return Continue(self._do_step)

    def finalize(self) -> t.Optional[ExitCode]:
        """Finalize the workgraph.
        Output the results of the workgraph and the new data.
        """
        # expose outputs of the workgraph
        self.task_manager.state_manager.update_meta_tasks("graph_ctx")
        self.wg.update()
        graph_results = {}
        # print("graph results:", self.ctx._task_results.get("_flow", {}))
        for name, data in self.ctx._task_results.get("_flow", {}).items():
            socket = self.wg.tasks[data["task_name"]].outputs[data["socket_name"]]
            graph_results[name] = (
                socket._value
                if isinstance(socket, TaskSocketNamespace)
                else socket.value
            )
        self.out_many(graph_results)
        # output the new data
        if self.ctx._new_data:
            self.out("new_data", self.ctx._new_data)
        self.report("Finalize workgraph.")
        for task in self.wg.tasks:
            if (
                self.task_manager.state_manager.get_task_runtime_info(
                    task.name, "state"
                )
                == "FAILED"
            ):
                return self.exit_codes.TASK_FAILED


# ---------------------------------------------------------------------------------------


async def wait_for(
    socket: TaskSocket, interval: float = 5.0, timeout: float = 604800.0
) -> None:
    """
    Wait for the socket's node to reach a terminal state, with a timeout.

    :param socket: The TaskSocket instance to monitor.
    :param interval: How often to check the node state.
    :param timeout: Maximum time to wait before raising a TimeoutError.
    """
    import time
    import asyncio

    start_time = time.monotonic()
    while socket._node.state not in ["FINISHED", "FAILED"]:
        if time.monotonic() - start_time > timeout:
            print("Timeout reached while waiting for node state.")
            return
        print(f"Task {socket._node.name} state:", socket._node.state)
        await asyncio.sleep(interval)
    socket._node.graph.update()
