"""
======================================
Recursive workflows with Graph Builder
======================================
"""
# %%
# For this how-to please familiarize yourself before with the concepts:
# - `Graph Builder <graph_builder.html>`_
# - `Context <../context.html>`_

# %%
# Introduction
# ============
# In this how-to, you will learn how to create a recursive workflow using
# the Graph Builder. As recursive functions invoke themself we need to define
# the recursive function in a separate module to be able execute it.
# TODO Xing maybe you can add some technical aspects. For me it is a bit unclear

# Load AiiDA profile
from aiida import load_profile

load_profile()


# %%
# Factorial example
# =================
# In this example we will define a factorial function using the Graph Builder
#
# .. code-block:: python
#
#    def factorial(n):
#        if n == 1:
#            return n
#        else:
#            return n * factorial(n-1)

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


# %%
# However we this alone will not work when we execute we will get an error that
# factorial is not defined.

try:
    wg = WorkGraph()
    wg.add_task(factorial, n=2)
    wg.run()
except Exception as err:
    print(err)

# %%
# We have to load the definition of the graph builder from a module. We copied the
# definiton above and put it into this file `factorial_graph_builder.py <factorial_graph_builder.py>`_.

from factorial_graph_builder import factorial

wg = WorkGraph()
wg.add_task(factorial, n=2)
wg.to_html()

# %%
# We run the workgraph
wg.submit(wait=True)

# %%
# We can also see the recursive execution of the workgraphs in the provenance graph

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Note that by each recursion steep you create a new WorkGraph and the limit of
# the number of processes (WorkGraph, CalcJobs, etc.) is 200 processes in AiiDA.
# It can as increased as explaind in
# `parallel how-to <parallel.html#maximum-number-of-active-workgraphs>.
