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

# %%
# Use ``aiida-core`` components in a ``WorkGraph``
# ------------------------------------------------
#
# ``aiida-core`` components are fully compatible with ``WorkGraph``.
# This means that any ``aiida-core`` component can be cast as a task and used within a ``WorkGraph``.
#
# In the following example, we combine all four ``aiida-core`` processes in a single ``WorkGraph``, connecting their inputs and outputs as needed.

# %%
# We start by loading the AiiDA profile and importing the necessary components:

from aiida import load_profile, orm

load_profile()

# %%
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida.engine import calcfunction, workfunction
from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

from aiida_workgraph import task, spec

# %%
# Next, let's define a ``calcfunction`` and ``workfunction``


@calcfunction
def add(x, y):
    return x + y


@workfunction
def add_more(x, y, z):
    sum_x_y = add(x, y)
    return add(sum_x_y, z)


# %%
# To use ``aiida-core`` components in a ``WorkGraph``, we can simply cast them as tasks using the ``task`` decorator functionally (``task(<aiida-core-component>)``).
# Task functional components can then be called functionally with their respective inputs, linking outputs to inputs as needed.


@task.graph
def AiiDAComponentsWorkflow():
    calcjob_sum = task(ArithmeticAddCalculation)(
        code=orm.load_code("add@localhost"),
        x=1,
        y=2,
    ).sum  # `ArithmeticAddCalculation` explicitly defines a 'sum' output

    calcfunction_sum = task(add)(
        x=calcjob_sum,
        y=4,
    ).result  # `calcfunction` implicitly defines a 'result' output

    workchain_result = task(MultiplyAddWorkChain)(
        code=orm.load_code("add@localhost"),
        x=calcfunction_sum,
        y=2,
        z=3,
    ).result  # `MultiplyAddWorkChain` explicitly defines a 'result' output

    workfunction_sum = task(add_more)(
        x=workchain_result,
        y=2,
        z=3,
    ).result  # `workfunction` implicitly defines a 'result' output

    return workfunction_sum


wg = AiiDAComponentsWorkflow.build_graph()
wg.to_html()

# %%
# .. important::
#
#    The ``task`` decorator accepts the same arguments when used functionally.
#    However, here are a few things to consider:
#
#    - For ``WorkChain`` and ``CalcJob``, inputs and outputs are determined by the AiiDA process definition and cannot be overridden.
#    - For ``calcfunction`` and ``workfunction``:
#
#      - When the AiiDA process provides a single output (e.g., ``return x``), the output socket name is implicitly set to ``result``.
#        The user may override this by assigning a custom output socket name (e.g., ``task(outputs=["my_sum"])(add)(...)``).
#        However, the benefits of clear graph visualization labels must be weighed against the loss of provenance consistency (the provenance graph will still show ``result``).
#      - When the AiiDA process provides multiple outputs (e.g., ``return {"x": x, "y": y}``), it is actually necessary to assign output socket names explicitly.
#        However, the previous point applies here as well if the user chooses to assign output socket names different from the dictionary keys of the AiiDA process.
#
#    For more about sockets, please refer to the :doc:`../../concept/autogen/socket_concept` concept section.


# %%
# Let's run our ``aiida-core``-powered ``WorkGraph`` and examine the provenance graph:

wg.run()

# %%
from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

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
        spec.input_namespace("workgraph", dynamic=True)
        spec.outline(
            cls.run_workgraph,
            cls.results,
        )
        spec.output("sum")
        spec.output("product")

    def run_workgraph(self):
        from aiida_workgraph.engine.workgraph import WorkGraphEngine

        process = self.submit(WorkGraphEngine, **self.inputs.workgraph)
        self.to_context(workgraph_process=process)

    def results(self):
        self.out("sum", self.ctx.workgraph_process.outputs.sum)
        self.out("product", self.ctx.workgraph_process.outputs.product)


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
def IntegratedAddMultiply() -> spec.namespace(sum=any, product=any):
    the_sum = add(1, 2).result
    the_product = multiply(the_sum, 3).result
    return {"sum": the_sum, "product": the_product}


wg = IntegratedAddMultiply.build_graph()

# %%
# We can export our workgraph as a dictionary using the ``prepare_inputs()`` method and use it as the input to our ``WorkChain``:

from aiida.engine import run_get_node

inputs = {"workgraph": wg.prepare_inputs(metadata={"call_link_label": "workgraph"})}
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
#       builder.workgraph_data = wg.prepare_inputs(metadata={"call_link_label": "workgraph"})
#
# Let's check the outputs of our ``WorkChain``:


print("Results:")
print("  Sum:    ", result["sum"])
print("  Product:", result["product"])

# %%
# And finally, we can have a look at the provenance graph:

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(node.pk)

# %%
# Further reading
# ---------------
# One can also use ``WorkGraph`` inside a ``WorkChain``, please refer to the `Calling WorkGraph within a WorkChain <workchain_call_workgraph.ipynb>`_ for more details.
