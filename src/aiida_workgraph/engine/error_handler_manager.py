from __future__ import annotations
import traceback


class ErrorHandlerManager:
    def __init__(self, process, ctx_manager, logger):
        self.process = process
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger

    def run_error_handlers(self, task_name: str) -> None:
        """Run error handlers for a task."""

        self.process.report(f"Run error handlers for {task_name}")

        node = self.process.task_manager.state_manager.get_task_runtime_info(
            task_name, "process"
        )
        if not node or not node.exit_status:
            return
        # error_handlers from the task
        error_handlers = self.process.node.task_error_handlers.get(task_name, {})
        for data in error_handlers.values():
            if node.exit_status in data.get("exit_codes", []):
                handler = data["handler"]
                self.run_error_handler(handler, data, task_name)
                return
        # error_handlers from the workgraph
        for data in self.process.wg._error_handlers.values():
            if node.exit_code.status in data["tasks"].get(task_name, {}).get(
                "exit_codes", []
            ):
                handler = data["handler"]
                metadata = data["tasks"][task_name]
                self.run_error_handler(handler, metadata, task_name)
                return

    def run_error_handler(self, handler: dict, metadata: dict, task_name: str) -> None:
        """Run the error handler for a task."""
        from inspect import signature
        from node_graph.executor import NodeExecutor

        handler = NodeExecutor(**handler).executor
        handler_sig = signature(handler)
        metadata.setdefault("retry", 0)
        self.process.report(f"Run error handler: {handler.__name__}")
        if metadata["retry"] < metadata["max_retries"]:
            task = self.process.task_manager.get_task(task_name)
            try:
                # Run the error handler to update the inputs of the task
                if "engine" in handler_sig.parameters:
                    msg = handler(task, engine=self, **metadata.get("kwargs", {}))
                else:
                    msg = handler(task, **metadata.get("kwargs", {}))
                # Reset the task to rerun it
                self.process.task_manager.state_manager.reset_task(task.name)
                if msg:
                    self.process.report(msg)
                metadata["retry"] += 1
            except Exception as e:
                error_traceback = traceback.format_exc()  # Capture the full traceback
                self.logger.error(
                    f"Error in running error handler for {task_name}: {e}\n{error_traceback}"
                )
                self.process.report(
                    f"Error in running error handler for {task_name}: {e}\n{error_traceback}"
                )
