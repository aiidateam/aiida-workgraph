from typing import Any


def get_context(context: dict, key: str) -> Any:
    """Get the context value."""
    results = {"result": context._task_results["graph_ctx"].get(key)}
    return results


def update_ctx(context: dict, key: str, value: Any) -> None:
    """Set the context value."""
    context._task_results["graph_ctx"][key] = value


def select(condition, true=None, false=None):
    """Select the data based on the condition."""
    if condition:
        return true
    return false
