import asyncio
from aiida_workgraph.task import SpecTask
from typing import Callable, Optional, Any, Dict
from node_graph.socket_spec import SocketSpec
from node_graph.node_spec import NodeSpec
from .function_task import build_callable_nodespec
from aiida_workgraph.socket_spec import namespace
from node_graph.executor import RuntimeExecutor
from node_graph.error_handler import ErrorHandlerSpec


class AwaitableFunctionTask(SpecTask):
    """Awaitable task with function as executor."""

    identifier = "workgraph.awaitable_function"
    name = "awaitable"
    node_type = "awaitable"
    catalog = "Control"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable
        if var_kwargs is None:
            awaitable_target = asyncio.ensure_future(
                executor(*args, **kwargs),
                loop=engine_process.loop,
            )
        else:
            awaitable_target = asyncio.ensure_future(
                executor(*args, **kwargs, **var_kwargs),
                loop=engine_process.loop,
            )
        return awaitable_target, "FINISHED"


class MonitorFunctionTask(SpecTask):
    """Monitor task with function as executor."""

    identifier = "workgraph.monitor_function"
    name = "monitor"
    node_type = "monitor"
    catalog = "Control"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.tasks.monitors import monitor

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable
        # get the raw function without the decorator
        if hasattr(executor, "_func"):
            executor = executor._func
        # add function and interval to the args
        args = [
            executor,
            kwargs.pop("interval", 1),
            kwargs.pop("timeout", 3600),
            *args,
        ]

        if var_kwargs is None:
            awaitable_target = asyncio.ensure_future(
                monitor(*args, **kwargs),
                loop=engine_process.loop,
            )
        else:
            awaitable_target = asyncio.ensure_future(
                monitor(*args, **kwargs, **var_kwargs),
                loop=engine_process.loop,
            )
        return awaitable_target, "FINISHED"


def _build_awaitable_function_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> NodeSpec:

    return build_callable_nodespec(
        obj=obj,
        node_type="AWAITABLE",
        base_class=AwaitableFunctionTask,
        identifier=identifier,
        process_cls=None,
        in_spec=in_spec,
        out_spec=out_spec,
        error_handlers=error_handlers,
    )


def _build_monitor_function_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> NodeSpec:
    # defaults for interval/timeout — set on the NAMESPACE (keys = field names)
    add_in = namespace(
        interval=(int, 5),
        timeout=(int, 3600),
    )
    add_out = namespace(exit_code=Any)

    return build_callable_nodespec(
        obj=obj,
        node_type="MONITOR",
        base_class=MonitorFunctionTask,
        identifier=identifier,
        process_cls=None,
        in_spec=in_spec,
        out_spec=out_spec,
        add_inputs=add_in,
        add_outputs=add_out,
        error_handlers=error_handlers,
    )
