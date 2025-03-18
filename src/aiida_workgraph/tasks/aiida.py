from aiida_workgraph.task import Task


class AiiDAFunctionTask(Task):
    """Task with AiiDA calcfunction/workfunction as executor."""

    identifier = "workgraph.aiida_functions"
    name = "aiida_function"
    node_type = "function"
    catalog = "AIIDA"

    def execute(self, args=None, kwargs=None, var_kwargs=None):
        from aiida.engine import run_get_node
        from node_graph.executor import NodeExecutor

        executor = NodeExecutor(**self.get_executor()).executor
        # the imported executor could be a wrapped function
        if hasattr(executor, "_NodeCls") and hasattr(executor, "_func"):
            executor = getattr(executor, "_func")
        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": self.name})
        # since aiida 2.5.0, we need to use args_dict to pass the args to the run_get_node
        if var_kwargs is None:
            _, process = run_get_node(executor, **kwargs)
        else:
            _, process = run_get_node(executor, **kwargs, **var_kwargs)
        process.label = self.name

        return process, "FINISHED"


class CalcFunctionTask(AiiDAFunctionTask):
    identifier = "workgraph.calcfunction"
    name = "calcfunction"
    node_type = "CalcFunction"
    catalog = "AIIDA"


class WorkFunctionTask(AiiDAFunctionTask):
    identifier = "workgraph.workfunction"
    name = "workfunction"
    node_type = "WorkFunction"
    catalog = "AIIDA"


class AiiDAProcessTask(Task):
    """Task with AiiDA calcfunction/workfunction as executor."""

    identifier = "workgraph.aiida_process"
    name = "aiida_process"
    node_type = "Process"
    catalog = "AIIDA"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from node_graph.executor import NodeExecutor
        from aiida_workgraph.utils import create_and_pause_process

        executor = NodeExecutor(**self.get_executor()).executor

        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": self.name})
        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                executor,
                kwargs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(executor, **kwargs)
            state = "RUNNING"
        process.label = self.name

        return process, state


class CalcJobTask(AiiDAProcessTask):
    identifier = "workgraph.calcjob"
    name = "calcjob"
    node_type = "CalcJob"
    catalog = "AIIDA"


class WorkChainTask(AiiDAProcessTask):
    identifier = "workgraph.workchain"
    name = "workchain"
    node_type = "WorkChain"
    catalog = "AIIDA"
