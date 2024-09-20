from aiida_workgraph import task, WorkGraph


@task.calcfunction
def identity(x):
    return x.clone()


@task.calcfunction
def multiply(x, y):
    return x * y


@task.graph_builder(outputs=[{"name": "result", "from": "context.task_out"}])
def factorial(n):
    wg = WorkGraph()
    if n == 1:
        task = wg.add_task(identity, name="identity", x=n)
    else:
        factorial_recursion = wg.add_task(factorial, name="factorial", n=n - 1)
        task = wg.add_task(
            multiply, name="multiply", x=n, y=factorial_recursion.outputs["result"]
        )

    # set the identity or multiply task as output
    task.set_context({"result": "task_out"})
    return wg
