from __future__ import annotations


class ErrorHandlerManager:
    def __init__(self, process, ctx_manager, logger):
        self.process = process
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger

    def run_error_handlers(self, task_name: str) -> None:
        """Run error handlers for a task."""

        node = self.process.task_manager.get_task_state_info(task_name, "process")
        if not node or not node.exit_status:
            return
        # error_handlers from the task
        for _, data in self.ctx._tasks[task_name]["error_handlers"].items():
            if node.exit_status in data.get("exit_codes", []):
                handler = data["handler"]
                self.run_error_handler(handler, data, task_name)
                return
        # error_handlers from the workgraph
        for _, data in self.ctx._error_handlers.items():
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
        from aiida_workgraph.utils import get_executor

        handler, _ = get_executor(handler)
        handler_sig = signature(handler)
        metadata.setdefault("retry", 0)
        self.process.report(f"Run error handler: {handler.__name__}")
        if metadata["retry"] < metadata["max_retries"]:
            task = self.process.task_manager.get_task(task_name)
            try:
                if "engine" in handler_sig.parameters:
                    msg = handler(task, engine=self, **metadata.get("kwargs", {}))
                else:
                    msg = handler(task, **metadata.get("kwargs", {}))
                self.process.task_manager.update_task(task)
                if msg:
                    self.process.report(msg)
                metadata["retry"] += 1
            except Exception as e:
                self.process.report(f"Error in running error handler: {e}")
