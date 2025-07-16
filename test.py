from aiida_workgraph import WorkGraph, task, If, While
from aiida_workgraph.utils import generate_node_graph
from aiida import load_profile

load_profile()

@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y

@task.graph_builder(outputs=["result"])
def add_multiply_if(x, y):
    # TODO: Context manager here
    wg = WorkGraph()
    if x.value < 0:
        add1 = wg.add_task(add, name="add1", x=x, y=y)
        # export the result of add1 to the graph-level outputs
        wg.outputs.result = add1.outputs.result
    else:
        multiply1 = wg.add_task(multiply, name="multiply1", x=x, y=y)
        # export the result of multiply1 to the graph-level outputs
        wg.outputs.result = multiply1.outputs.result
    return wg

# TODO: Context manager here
with WorkGraph("if_graph_builer") is wg:
    sum = add(x=1, y=1)
    add_multiply_if1 = wg.add_task(
        add_multiply_if, name="add_multiply_if1", x=add1.outputs.result, y=3
    )
    add1 = wg.add_task(add, name="add2", x=add_multiply_if1.outputs.result, y=1)

import ipdb; ipdb.set_trace()

wg.run()
