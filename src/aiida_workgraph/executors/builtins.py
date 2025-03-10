from typing import Any


def get_context(context: dict, key: str) -> Any:
    """Get the context value."""
    results = {"result": getattr(context, key)}
    return results


def set_context(context: dict, key: str, value: Any) -> None:
    """Set the context value."""
    setattr(context, key, value)


def select(condition, true=None, false=None):
    """Select the data based on the condition."""
    if condition:
        return true
    return false
