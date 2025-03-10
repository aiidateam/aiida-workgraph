import asyncio
from aiida_workgraph import Task
from .function_task import DecoratedFunctionTaskFactory


class AwaitableFunctionTask(Task):
    """Awaitable task with function as executor."""

    identifier = "workgraph.awaitable_function"
    name = "awaitable"
    node_type = "awaitable"
    catalog = "Control"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from node_graph.executor import NodeExecutor

        executor = NodeExecutor(**self.get_executor()).executor
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


class MonitorFunctionTask(Task):
    """Monitor task with function as executor."""

    identifier = "workgraph.monitor_function"
    name = "monitor"
    node_type = "monitor"
    catalog = "Control"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from node_graph.executor import NodeExecutor
        from aiida_workgraph.executors.monitors import monitor

        executor = NodeExecutor(**self.get_executor()).executor
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


class AwaitableFunctionTaskFactory(DecoratedFunctionTaskFactory):
    """A factory to create an AwaitableFunctionTask from functions."""

    default_task_type = "awaitable"
    default_base_class = AwaitableFunctionTask


class MonitorFunctionTaskFactory(DecoratedFunctionTaskFactory):
    """A factory to create a MonitorFunctionTask from functions."""

    default_task_type = "monitor"
    default_base_class = MonitorFunctionTask
