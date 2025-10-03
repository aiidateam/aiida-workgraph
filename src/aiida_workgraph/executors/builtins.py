from typing import Any


def UnavailableExecutor(*args, **kwargs):
    raise RuntimeError('This executor was defined dynamically and is not available from the database snapshot.')


def get_context(context: dict, key: str) -> Any:
    """Get the context value."""
    results = {'result': context._task_results['graph_ctx'].get(key)}
    return results


def update_ctx(context: dict, key: str, value: Any) -> None:
    """Set the context value."""
    context._task_results['graph_ctx'][key] = value


def select(condition, true=None, false=None):
    """Select the data based on the condition."""
    if condition:
        return true
    return false


def get_item(data: dict, key: str) -> Any:
    """Get an item from a dictionary."""
    return data.get(key, None)


def return_input(**kwargs: Any) -> dict:
    """Return the input"""
    return kwargs
