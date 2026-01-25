"""
Interoperate with ``aiida-core`` components
===========================================
"""

# %%
# Introduction
# ------------
#
# If you’re already familiar with AiiDA, you may be interested in integrating ``aiida-core`` components, such as ``CalcJob``, ``calcfunction``, ``WorkChain``, ``workfunction``, and ``ProcessBuilder`` within a ``WorkGraph``.
# This integration enables you to seamlessly reuse existing workflows built with these paradigms inside a more flexible graph-based structure.
# Similarly, if you have existing robust workflows built with ``aiida-core`` components, you can easily incorporate ``WorkGraph`` components into them.
#
# This guide will demonstrate the two-way integration of ``WorkGraph`` with ``aiida-core`` components, providing examples for each.
#
# .. note::
#
#    This guide assumes prior knowledge of ``aiida-core`` components.
#    If you’re unfamiliar with them, please refer to the official documentation on `Calculations <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/index.html>`_ and `Workflows <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/workflows/index.html>`_.

import typing as t
from aiida import load_profile, orm
from aiida.engine import calcfunction
from aiida_workgraph import namespace, task

load_profile()

# %%
# Use ``aiida-core`` components in a ``WorkGraph``
# ------------------------------------------------
#
# ``aiida-core`` components are fully compatible with ``WorkGraph``.
# This means that any ``aiida-core`` component can be cast as a task and used within a ``WorkGraph``.
#
# In the following example, we combine all four ``aiida-core`` processes in a single ``WorkGraph``, connecting their inputs and outputs as needed.
#
# Calcfunction and workfunction
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Let's first define a ``calcfunction`` and ``workfunction``
from aiida import orm


@calcfunction
def add_multiply(x, y):
    return {'sum': orm.Int(x + y), 'product': orm.Int(x * y)}


# %%
# To use ``aiida-core`` components in a ``WorkGraph``, we can simply cast them as tasks using the ``task`` decorator functionally (``task(<aiida-core-component>)``).

add_multiply_task = task(outputs=namespace(sum=int, product=int))(add_multiply)

# %%
# One can also use the decorator syntax directly


@task.calcfunction(outputs=namespace(sum=int, product=int))
def add_multiply_task(x, y):
    return {'sum': orm.Int(x + y), 'product': orm.Int(x * y)}


# %%
# The same syntax applies to `workfunction`.


# %%
#
# .. important::
#
#    The ``task`` decorator accepts the same arguments when used functionally.
#    However, here are a few things to consider:
#      - Input specifications are inferred from the Python signature. Declaring explicit input namespaces is not supported, because `calcfunction` and `workfunction` accept namespaced inputs only via keyword argument (e.g., ``**kwargs``): this is treated as a dynamic namespace and bypasses validation, allowing arbitrary nested AiiDA data to be provided.
#      - Outputs are always dynamic: whatever structure the function returns becomes its provenance without validation. Defining an output specification is therefore only needed to expose output sockets that other tasks can link to.
#      - When the AiiDA process provides a single output (e.g., ``return x``), the socket name defaults to ``result``. You may override it (e.g., ``task(outputs=["my_sum"])(add)(...)``) for clearer graph labels, understanding the provenance graph will still record ``result``.
#      - When the process returns a dictionary (e.g., ``return {"x": x, "y": y}``), you must assign output socket names explicitly if you want to reference individual entries. See :doc:`../../gallery/howto/autogen/annotate_inputs_outputs` for more on namespaces and dynamic sockets.
#
# CalcJob and WorkChain
# ^^^^^^^^^^^^^^^^^^^^^^^
# Similarly, ``CalcJob`` and ``WorkChain`` can also be cast as tasks using the ``task`` decorator.

from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation


ArithmeticAddTask = task(ArithmeticAddCalculation)
MultiplyAddTask = task(MultiplyAddWorkChain)


# %%
# .. important::
#
#    For ``WorkChain`` and ``CalcJob``, inputs and outputs are determined by the AiiDA process definition and cannot be overridden by passing the `inputs` or `outputs` arguments to the ``task`` decorator.
#

# %%
# Task components can then be called functionally with their respective inputs, linking outputs to inputs as needed.


def add(x, y):
    return x + y


def add_more(x, y, z):
    sum_x_y = add(x, y)
    return add(sum_x_y, z)


AddTask = task(add)
AddMoreTask = task(add_more)


@task.graph
def AiiDAComponentsWorkflow():
    calcjob_sum = ArithmeticAddTask(
        code=orm.load_code('add@localhost'),
        x=1,
        y=2,
    ).sum  # `ArithmeticAddCalculation` explicitly defines a 'sum' output

    calcfunction_sum = AddTask(
        x=calcjob_sum,
        y=4,
    ).result  # `calcfunction` implicitly defines a 'result' output

    workchain_result = MultiplyAddTask(
        code=orm.load_code('add@localhost'),
        x=calcfunction_sum,
        y=2,
        z=3,
    ).result  # `MultiplyAddWorkChain` explicitly defines a 'result' output

    workfunction_sum = AddMoreTask(
        x=workchain_result,
        y=2,
        z=3,
    ).result  # `workfunction` implicitly defines a 'result' output

    return workfunction_sum


wg = AiiDAComponentsWorkflow.build()
wg

# %%
# Let's run our ``aiida-core``-powered ``WorkGraph`` and examine the provenance graph:

wg.run()

# %%


wg.generate_provenance_graph()

# %%
# Use ``WorkGraph`` in ``WorkChain``
# ----------------------------------
#
# ``WorkGraph`` can also be used within ``aiida-core`` components.
# Whether you want to integrate a ``WorkGraph`` into an existing robust ``WorkChain``, or simply prefer to keep certain tasks as ``aiida-core`` components while using ``WorkGraph`` for others, incorporating ``WorkGraph`` into ``aiida-core`` components is a straightforward process.
#
# Let's define a ``WorkChain`` that submits a ``WorkGraph``:

from aiida.engine import WorkChain


class TestWorkChain(WorkChain):
    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input_namespace('workgraph', dynamic=True)
        spec.outline(
            cls.run_workgraph,
            cls.results,
        )
        spec.output('sum')
        spec.output('product')

    def run_workgraph(self):
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        process = self.submit(WorkGraphEngine, **self.inputs.workgraph)
        self.to_context(workgraph_process=process)

    def results(self):
        self.out('sum', self.ctx.workgraph_process.outputs.sum)
        self.out('product', self.ctx.workgraph_process.outputs.product)


# %%
# A few things to note about the above ``WorkChain``:
#
# - ``WorkGraphEngine`` is the AiiDA process associated with workgraphs.
#   Its ``workgraph_data`` expects a ``WorkGraph`` in dictionary format.
# - When gathering results, we expect to find workgraph outputs named ``sum`` and ``product``.
#
# Let's define a simple *AddMultiply* workgraph to provide as input to our ``WorkChain``:


@task
def add(x, y):
    return x + y


@task
def multiply(x, y):
    return x * y


@task.graph
def IntegratedAddMultiply() -> t.Annotated[dict, namespace(sum=int, product=int)]:
    the_sum = add(1, 2).result
    the_product = multiply(the_sum, 3).result
    return {'sum': the_sum, 'product': the_product}


wg = IntegratedAddMultiply.build()

# %%
# We can export our workgraph as a dictionary using the ``to_engine_inputs()`` method and use it as the input to our ``WorkChain``:

from aiida.engine import run_get_node

inputs = {'workgraph': wg.to_engine_inputs(metadata={'call_link_label': 'workgraph'})}
result, node = run_get_node(TestWorkChain, **inputs)

# %%
#
# .. tip::
#
#    All AiiDA processes classes, including ``WorkGraphEngine``, offer a ``get_builder()`` method.
#    You can use this method to extract the associated ``ProcessBuilder`` and use it to set inputs directly.
#
#    .. code:: python
#
#       builder = WorkGraphEngine.get_builder()
#       builder.workgraph_data = wg.to_engine_inputs(metadata={"call_link_label": "workgraph"})
#
# Let's check the outputs of our ``WorkChain``:


print('Results:')
print('  Sum:    ', result['sum'])
print('  Product:', result['product'])

# %%
# And finally, we can have a look at the provenance graph:
from aiida_workgraph.utils import generate_provenance_graph

generate_provenance_graph(node.pk)

# %%
# Further reading
# ---------------
# One can also use ``WorkGraph`` inside a ``WorkChain``, please refer to the `Calling WorkGraph within a WorkChain <workchain_call_workgraph.ipynb>`_ for more details.
