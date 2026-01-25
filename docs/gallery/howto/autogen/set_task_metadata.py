"""
Set task metadata
=================

By default, tasks are named after their function or process class name.
When a workflow becomes complex, or when the same task is used multiple times, giving tasks custom names and descriptions can greatly improve clarity.

This example shows how to set three useful metadata fields on tasks:

- ``call_link_label``: the **task name** in the WorkGraph, and the **link label** shown in the provenance visualization.
- ``label``: the **process label** stored on the AiiDA process.
- ``description``: a **human friendly description** stored on the AiiDA process.

We also show how to query the provenance by these fields.

.. note::

   For ``CalcJob`` and ``PythonJob`` tasks there are more metadata options.
   Please refer to the section on :ref:`PythonJob metadata <pythonjob_metadata>` and the `CalcJob options <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/usage.html#options>`_ for more details.

"""
# %%
# Here's an example:

from aiida import load_profile
from aiida_workgraph import task

load_profile()


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


@task.graph()
def add_and_multiply(x, y, z):
    sum_result = add(
        x,
        y,
        metadata={
            'call_link_label': 'my_add',
            'label': 'AdditionStep',
            'description': 'Compute x + y and store it as an intermediate result.',
        },
    ).result

    product_result = multiply(
        sum_result,
        z,
        metadata={
            'call_link_label': 'my_multiply',
            'label': 'MultiplicationStep',
            'description': 'Multiply the sum by z to produce the final result.',
        },
    ).result

    return product_result


# Build and visualize the graph
wg = add_and_multiply.build(x=1, y=2, z=3)
wg

# %%
# As you can see in the visualization, the task are now labeled "my_add" and "my_multiply" instead of the default function names.
#
# Run the workflow:

wg.run()


# %%
# Query provenance by label or description
# ----------------------------------------
# You can later find these processes by their metadata.
from aiida import orm

qb = orm.QueryBuilder()
qb.append(
    orm.ProcessNode,
    filters={'label': {'==': 'MultiplicationStep'}},
    project=['uuid', 'label', 'description'],
)
matches = qb.all()

print(f"Found {len(matches)} process node(s) labeled 'MultiplicationStep'")
for uuid, label, description in matches:
    print(f'- {uuid} | {label} | {description}')
