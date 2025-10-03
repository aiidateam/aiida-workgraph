from aiida_workgraph import task


@task
def add(x, y):
    """Add two numbers."""
    return x + y


@task
def multiply(x, y):
    """Multiply two numbers."""
    return x * y
