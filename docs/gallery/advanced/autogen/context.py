"""
=========================================
Use context variables in WorkGraphs
=========================================
"""

# %%
# In AiiDA workflows (both, traditional ``WorkChain`` s, as well as the ``WorkGraph``), the **Context** (typically
# represented by ``ctx``), is an internal container that can hold data shared between different tasks.
# It's particularly useful for more complex workflows.

# %%
# When to use
# -----------
#
# 1. Many-to-many relationships where direct linking becomes unwieldy
# 2. Conditional workflows where links depend on runtime conditions
# 3. Graph builders that need to expose selective internal state

# %%
# Setup
# -----

from aiida_workgraph import WorkGraph, task
from aiida import load_profile
from aiida.orm import Int

load_profile()

# %%
# Passing data to the context
# ---------------------------
#
# There are three ways to set data in the **Context**:
#
# - Setting the ``ctx`` attribute of the WorkGraph directly
# - Using the ``update_ctx`` method of the WorkGraph
# - Using the ``workgraph.set_context`` task
#
# Let's have a look how to use each of them:

wg1 = WorkGraph(name='context_1')

# 1. Setting the ``ctx`` attribute of the WorkGraph directly, on initialization
wg1.ctx = {'x': Int(2), 'data.y': Int(3)}

# 2. Using the ``update_ctx`` method
wg2 = WorkGraph(name='context_2')


@task.calcfunction()
def add(x, y):
    return x + y


add1 = wg2.add_task(add, 'add1', x=2, y=3)

# set result of add1 to ctx.sum
wg2.update_ctx({'sum': add1.outputs.result})

# 3. Using the ``workgraph.set_context`` task to set either a task result (socket) or a resolved value to the ctx
wg3 = WorkGraph(name='context_3')
add1 = wg3.add_task(add, 'add1', x=2, y=3)
wg3.add_task('workgraph.set_context', name='set_ctx1', key='sum', value=add1.outputs.result)

# %%
# Nested context keys
# -------------------
# To organize the context data in a hierarchical structure, the keys may contain dots ``.``` that create nesting
# Here is an example, to group the results of multipl add tasks to `ctx.sum`:

wg = WorkGraph(name='ctx_nested')
add1 = wg.add_task(add, 'add1', x=1, y=2)
add2 = wg.add_task(add, 'add2', x=3, y=4)

wg.update_ctx({'sum.add1': add1.outputs.result})
wg.update_ctx({'sum.add2': add2.outputs.result})

# Or, alternatively:
# wg.update_ctx({
#     "sum": {
#         "add1": add1.outputs.result,
#         "add2": add2.outputs.result
#     }
# })

print(wg.ctx.sum)

# %%
# Use data from the context
# -------------------------
#
# Also for accessing data from the context, there are different approaches:

# %%
#
# 1. One can use elements from the ``wg.ctx.x``, directly, e.g., as inputs for other tasks
#
add2 = wg1.add_task(add, 'add2', x=wg1.ctx.x, y=3)

# 2. Similarly, nested context keys can be accessed, such as ``wg2.ctx.sum.add1``
print(wg.ctx.sum.add1)
# Gives: SocketAny(name='add1', value=None)

# 3. One can use the `get_context` task to get the data from ctx. **This task will be shown in the GUI**
#
wg3.add_task('workgraph.get_context', name='get_ctx1', key='sum.add1')

wg.show()
wg.to_html()


@task.graph()
def internal_add(x, y):
    outputs1 = add(x=x, y=y)
    outputs2 = add(x=outputs1.sum, y=5)
    return outputs2.sum


# Usage in a main WorkGraph
wg_main = WorkGraph('main')
builder_task = wg_main.add_task(internal_add, x=10, y=20)
final_task = wg_main.add_task(add, x=builder_task.outputs.result, y=100)

# In this example, the context can serve as a bridge between the internal workings of a subgraph and its external
# interface. This allows selectively exposing internal results as clean, named outputs.
#

# %%
# Context Tasks in the GUI
# ------------------------
#
# This example shows how context operations appear in the workflow visualization
# when using explicit context tasks instead of direct context access.
#
# .. note::
#
#    Context tasks make data flow visible in the GUI, which is helpful for
#    debugging and understanding workflow execution. Use them when you want
#    explicit visibility of context operations.

# %%

wg = WorkGraph(name='context_gui_demo')

# Set initial context values
wg.ctx = {'x': 2, 'multiplier': 10}

# Use context tasks - these appear as nodes in the GUI
get_x = wg.add_task('workgraph.get_context', name='get_x', key='x')
get_multiplier = wg.add_task('workgraph.get_context', name='get_multiplier', key='multiplier')

# Perform calculation using context values
add1 = wg.add_task(add, 'add1', x=get_x.outputs.result, y=get_multiplier.outputs.result)

# Store result back to context - also appears in GUI
wg.add_task(
    'workgraph.set_context',
    name='store_result',
    key='final_result',
    value=add1.outputs.result,
)

wg.to_html()
wg.show()

# %%
# As shown in the GUI, the ``get_context`` and ``set_context`` tasks appear as
# visible nodes in the workflow graph. This makes the data flow through context
# explicit and traceable, unlike direct context access via ``wg.ctx.x`` which
# is invisible in the visualization.

# %%
wg.run()
print('State of WorkGraph         : {}'.format(wg.state))
print('Result of add1            : {}'.format(wg.tasks.add1.outputs.result.value))

# %%

# Generate node graph from the AiiDA process


wg.generate_provenance_graph()

# %%
# .. note::
#
#    If you pass data from one task to another task through the context, you may need to use ``wait`` to wait for
#    the data to be ready. See :doc:`/howto/autogen/control_task_execution_order`.
