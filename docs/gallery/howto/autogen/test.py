from aiida_workgraph import WorkGraph, task
from aiida import load_profile

load_profile()

@task.calcfunction
def add(x, y):
    return x + y


@task.calcfunction
def multiply(x, y):
    return x * y

@task.calcfunction
def generate_random_number():
    return 42
import ipdb; ipdb.set_trace()

with WorkGraph("AddMultiply") as add_multiply:
    add_multiply.add_task(add)
    add_multiply.add_task(multiply, x=add_multiply.tasks.add.outputs.result)

...

with WorkGraph("AddMultiplyComposed") as add_multiply_composed:
    add_multiply_composed.add_task(generate_random_number)
    add_multiply_composed.add_task(add_multiply, name=add_multiply.name)
    # add_multiply_composed.add_link(
    #     add_multiply_composed.tasks.generate_random_number.outputs.result,
    #     add_multiply_composed.tasks.AddMultiply.tasks.multiply.inputs.y,
    # )

with WorkGraph("AddMultiplyExtension") as add_multiply_extended:
    add_multiply_extended.add_task(generate_random_number)
    add_multiply_extended.extend(add_multiply)
    add_multiply_extended.add_link(
        add_multiply_extended.tasks.generate_random_number.outputs.result,
        add_multiply_extended.tasks.multiply.inputs.y,
    )