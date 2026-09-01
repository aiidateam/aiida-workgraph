from aiida_workgraph.task import Task
from aiida_workgraph.enums import TaskAction, TaskState
from aiida.engine import Process
from typing import Callable, Optional, Dict
from node_graph.socket_spec import SocketSpec
from node_graph.task_spec import TaskSpec, SchemaSource
from .function_task import build_callable_TaskSpec
from node_graph.executor import RuntimeExecutor
from node_graph.error_handler import ErrorHandlerSpec, normalize_error_handlers
from aiida_workgraph.utils import inspect_aiida_component_type
from aiida_workgraph.socket_spec import from_aiida_process


class AiiDAFunctionTask(Task):
    """Task with AiiDA calcfunction/workfunction as executor."""

    identifier = 'workgraph.aiida_functions'
    name = 'aiida_function'
    task_type = 'function'
    catalog = 'AIIDA'

    def execute(self, engine_process=None, args=None, kwargs=None, var_kwargs=None):
        from aiida.engine import run_get_node
        from node_graph.task_spec import BaseHandle

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable
        # the imported executor could be a wrapped function
        if isinstance(executor, BaseHandle) and hasattr(executor, '_callable'):
            executor = getattr(executor, '_callable')
        kwargs.setdefault('metadata', {})
        kwargs['metadata'].update({'call_link_label': self.name})
        # since aiida 2.5.0, we need to use args_dict to pass the args to the run_get_node
        try:
            if var_kwargs is None:
                _, process = run_get_node(executor, **kwargs)
            else:
                _, process = run_get_node(executor, **kwargs, **var_kwargs)
        except Exception:
            # `run_get_node` stores and CALL-links the process node before running it,
            # then re-raises the wrapped function's exception straight through instead
            # of returning the node, so the caller loses the reference on failure.
            # Recover it from provenance via the call link label set above. A task
            # inside a While zone reuses the same label on every iteration, so several
            # CALL links can share it; take the most recently stored node. Recovery is
            # best-effort: it must never replace the exception being handled, so any
            # failure here (including no matching link at all, e.g. under
            # metadata={'store_provenance': False}) falls through to the bare `raise`
            # and the original exception propagates untouched.
            process = None
            if engine_process is not None:
                try:
                    linked = engine_process.node.base.links.get_outgoing(link_label_filter=self.name).all()
                    if linked:
                        process = max((link.node for link in linked), key=lambda node: node.pk)
                except Exception:
                    process = None
            if process is None:
                raise
            return process, TaskState.FAILED

        return process, TaskState.FINISHED


class AiiDAProcessTask(Task):
    """Task with AiiDA calcfunction/workfunction as executor."""

    identifier = 'workgraph.aiida_process'
    name = 'aiida_process'
    task_type = 'Process'
    catalog = 'AIIDA'

    @classmethod
    def build(
        cls,
        callable,
        attached_error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    ):
        attached_error_handlers = normalize_error_handlers(attached_error_handlers)
        in_spec, out_spec = from_aiida_process(callable)
        return TaskSpec(
            identifier=callable.__name__,
            schema_source=SchemaSource.CALLABLE,
            catalog='AIIDA',
            inputs=in_spec,
            outputs=out_spec,
            executor=RuntimeExecutor.from_callable(callable),
            attached_error_handlers=attached_error_handlers,
            base_class=cls,
            task_type=inspect_aiida_component_type(callable),
        )

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable

        kwargs.setdefault('metadata', {})
        kwargs['metadata'].update({'call_link_label': self.name})
        if self.action == TaskAction.PAUSE:
            engine_process.report(f'Task {self.name} is created and paused.')
            process = create_and_pause_process(
                engine_process.runner,
                executor,
                kwargs,
                state_msg='Paused through WorkGraph',
            )
            state = TaskState.CREATED
            process = process.node
        else:
            process = engine_process.submit(executor, **kwargs)
            state = TaskState.RUNNING

        return process, state


class CalcJobTask(AiiDAProcessTask):
    identifier = 'workgraph.calcjob'
    name = 'calcjob'
    task_type = 'CalcJob'
    catalog = 'AIIDA'


class WorkChainTask(AiiDAProcessTask):
    identifier = 'workgraph.workchain'
    name = 'workchain'
    task_type = 'WorkChain'
    catalog = 'AIIDA'


def _build_aiida_function_taskspec(
    obj: Callable,
    identifier: Optional[str] = None,
    catalog: str = 'AIIDA',
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> TaskSpec:
    from aiida_workgraph.utils import inspect_aiida_component_type
    from dataclasses import replace

    spec = build_callable_TaskSpec(
        obj=obj,
        task_type=inspect_aiida_component_type(obj),
        catalog=catalog,
        base_class=AiiDAFunctionTask,
        identifier=identifier,
        process_cls=Process,
        in_spec=in_spec,
        out_spec=out_spec,
        error_handlers=error_handlers,
    )
    # the outputs of calcfunctions/workfunctions are always dynamic
    spec = replace(spec, outputs=replace(spec.outputs, meta=replace(spec.outputs.meta, dynamic=True)))
    if obj.spec().inputs.dynamic:
        spec = replace(spec, inputs=replace(spec.inputs, meta=replace(spec.inputs.meta, dynamic=True)))
    return spec
