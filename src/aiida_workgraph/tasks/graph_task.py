from aiida_workgraph.task import Task
from typing import Callable, Optional
from node_graph.socket_spec import SocketSpec
from node_graph.node_spec import NodeSpec
from .function_task import build_callable_nodespec
from node_graph.executor import RuntimeExecutor
from aiida.engine import Process


class GraphTask(Task):
    """Graph builder task"""

    identifier = 'workgraph.graph_task'
    name = 'graph_task'
    node_type = 'graph_task'
    catalog = 'builtins'

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process, call_depth_from_node
        from aiida_workgraph.engine.workgraph import WorkGraphEngine
        from aiida_workgraph import task, WorkGraph
        from node_graph.utils.graph import materialize_graph
        from aiida_workgraph.task import TaskHandle

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable
        max_depth = self.spec.metadata.get('max_depth', 100)
        metadata = kwargs.pop('metadata', {}) if kwargs else {}
        metadata.setdefault('call_link_label', self.name)
        # Cloudpickle doesn’t restore the function’s own name in its globals after unpickling,
        # so any recursive calls would raise NameError. As a temporary workaround, we re-insert
        # the decorated function into its globals under its original name.
        # Downside: this mutates the module globals at runtime, if another symbol with the same name exists,
        # we may introduce hard-to-trace bugs or collisions.
        if isinstance(executor, TaskHandle) and hasattr(executor, '_callable'):
            executor = executor._callable
        # We override the executor's globals to include the decorated function for recursion
        executor.__globals__[executor.__name__] = task.graph(max_depth=max_depth)(executor)
        depth = call_depth_from_node(engine_process.node)
        if depth >= max_depth:
            if depth >= max_depth:
                msg = (
                    f"Graph task '{self.name}' exceeded the recursion safeguard.\n"
                    f'- Current AiiDA process call depth (approx.): {depth}\n'
                    f'- Allowed maximum          :                  {max_depth}\n'
                    f'- Process UUID:                               {engine_process.node.uuid}\n\n'
                    f'Deeply nested process calls (>100) are generally discouraged. '
                    f'Prefer wrapping iterative logic inside a single task instead of '
                    f'recursively spawning new graph tasks.\n\n'
                    f'However, if you are confident that recursion is the right design, '
                    f'you can explicitly set a higher limit in your decorator, e.g.:\n'
                    f'    @task.graph(max_depth=200)\n'
                )
                engine_process.report(msg)
                raise RecursionError(msg)
        wg = materialize_graph(
            executor,
            self.spec.inputs,
            self.spec.outputs,
            self.name,
            WorkGraph,
            args=args,
            kwargs=kwargs,
            var_kwargs=var_kwargs,
        )
        wg.parent_uuid = engine_process.node.uuid
        inputs = wg.to_engine_inputs(metadata=metadata)
        if self.action == 'PAUSE':
            engine_process.report(f'Task {self.name} is created and paused.')
            process = create_and_pause_process(
                engine_process.runner,
                WorkGraphEngine,
                inputs,
                state_msg='Paused through WorkGraph',
            )
            state = 'CREATED'
            process = process.node
        else:
            process = engine_process.submit(WorkGraphEngine, **inputs)
            state = 'RUNNING'

        return process, state


def _build_graph_task_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
    max_depth: int = 100,
    catalog: str = 'Others',
) -> NodeSpec:
    # defaults for max depth
    metadata = {'max_depth': max_depth}
    # We use Process as the process class here, so that the task inherits the metadata
    # inputs from the base Process class, such as 'call_link_label'.
    # While the actual process class will be the WorkGraphEngine,
    # which is set at runtime in the execute() method

    return build_callable_nodespec(
        obj=obj,
        node_type='GRAPH',
        base_class=GraphTask,
        identifier=identifier,
        catalog=catalog,
        process_cls=Process,
        in_spec=in_spec,
        out_spec=out_spec,
        metadata=metadata,
    )
