from aiida_workgraph.task import SpecTask
from typing import Callable, Optional
from node_graph.socket_spec import SocketSpec
from node_graph.node_spec import NodeSpec
from .function_task import build_callable_nodespec


class GraphTask(SpecTask):
    """Graph builder task"""

    identifier = "workgraph.graph_task"
    name = "graph_task"
    node_type = "graph_task"
    catalog = "builtins"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_workgraph.engine.workgraph import WorkGraphEngine
        from aiida_workgraph import task, WorkGraph
        from node_graph.utils.graph import materialize_graph
        from node_graph.node_spec import BaseHandle

        executor = self.get_executor().callable
        # Cloudpickle doesn’t restore the function’s own name in its globals after unpickling,
        # so any recursive calls would raise NameError. As a temporary workaround, we re-insert
        # the decorated function into its globals under its original name.
        # Downside: this mutates the module globals at runtime, if another symbol with the same name exists,
        # we may introduce hard-to-trace bugs or collisions.
        if isinstance(executor, BaseHandle) and hasattr(executor, "_func"):
            executor = executor._func
        if executor.__name__ not in executor.__globals__:
            executor.__globals__[executor.__name__] = task.graph()(executor)
        wg = materialize_graph(
            executor,
            self._spec.inputs,
            self._spec.outputs,
            self.name,
            WorkGraph,
            args=args,
            kwargs=kwargs,
            var_kwargs=var_kwargs,
        )
        wg.parent_uuid = engine_process.node.uuid
        inputs = wg.prepare_inputs(metadata={"call_link_label": self.name})
        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                WorkGraphEngine,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(WorkGraphEngine, **inputs)
            state = "RUNNING"
        process.label = self.name

        return process, state


def _build_graph_task_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
) -> NodeSpec:

    return build_callable_nodespec(
        obj=obj,
        node_type="GRAPH",
        base_class=GraphTask,
        identifier=identifier,
        process_cls=None,
        in_spec=in_spec,
        out_spec=out_spec,
    )
