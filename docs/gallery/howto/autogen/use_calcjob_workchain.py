"""
=======================================
Interoperate with aiida-core components
=======================================
"""

# %%
# Introduction
# ============
# If you’re already familiar with AiiDA, you may be interested in how to integrate ``aiida-core`` components—such as ``CalcJob``, ``calcfunction``, ``WorkChain``, ``workfunction``, and ``ProcessBuilder`` within a ``WorkGraph``. This integration enables you to seamlessly reuse existing workflows built with these paradigms inside a more flexible graph-based structure.
#
# Adding these components as tasks in an aiida-workgraph is straightforward, and their inputs and outputs can be easily connected to form complex workflows.
#
# .. note::
#
#    This guide assumes prior knowledge of ``aiida-core`` components. If you’re unfamiliar with them, refer to the official documentation on `Calculations <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/index.html>`_ and `Workflows <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/workflows/index.html>`_.

# %%
# Setting up AiiDA environment 
# ----------------------------
# We need to initialize the computer and code resources to establish where the workflow is run.

from aiida import load_profile

load_profile()

from aiida import orm
from aiida.common.exceptions import NotExistent

try:
    bash_code = orm.load_code(
        "bash@localhost"
    )  # The computer label can also be omitted here
except NotExistent:
    bash_code = orm.InstalledCode(
        label="bash",
        computer=orm.load_computer("localhost"),
        filepath_executable="/bin/bash",
        default_calc_job_plugin="core.arithmetic.add",
    ).store()

# %%
# Integrating ``aiida-core`` components
# =====================================
# We begin by examining how the ``aiida-core`` components are visualized in the task view. The following examples show that key information—such as inputs and outputs—is clearly displayed for each component.


from aiida_workgraph import WorkGraph

wg = WorkGraph("aiida_components")

# %%
# CalcJob
# -------
# As ``CalcJob`` we use the ``ArithmeticAddCalculation`` that adds two numbers ``x+y``.

from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

task_calcjob = wg.add_task(ArithmeticAddCalculation)
task_calcjob.to_html()

# %%
# WorkChain
# ---------
# As ``WorkChain`` we use the ``MultiplyAddWorkChain`` that performs the operation ``x*y+z``.

from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain

task_workchain = wg.add_task(MultiplyAddWorkChain)
task_workchain.to_html()

# %%
# calcfunction
# ------------
# We create a ``calcfunction`` that adds two numbers ``x+y``.

from aiida.engine import calcfunction


@calcfunction
def add(x, y):
    return x + y


task_workfunction = wg.add_task(add)
task_workfunction.to_html()

# %%
# workfunction
# ------------
# We create a ``workfunction`` that adds three numbers ``x+y+z``.

from aiida.engine import workfunction


@workfunction
def add_more(x, y, z):
    sum_x_y = add(x, y)
    return add(sum_x_y, z)


task_workfunction = wg.add_task(add_more, name="workfunction")
task_workfunction.to_html()


# %%
# Putting all together
# --------------------
# Now we create a workflow using all of these components, linking their inputs and outputs just as in earlier ``WorkGraph`` examples.

wg = WorkGraph("integrate_aiida_core_components")

# We demonstrate how to use the process builder to set inputs for the first task
builder = ArithmeticAddCalculation.get_builder()

# Assign code and inputs to the builder
builder.code = bash_code
builder.x = orm.Int(1)
builder.y = orm.Int(2)

# 1+2 = 3
add_calcjob_task = wg.add_task(ArithmeticAddCalculation, name="add_calcjob")
add_calcjob_task.set_from_builder(builder)

# Alternatively, we can initalize the ``add_calcjob`` task without the process builder
# add_calcjob_task = wg.add_task(
#     ArithmeticAddCalculation, name="add_calcjob", x=1, y=2, code=bash_code
# )

# 3+4 = 7
wg.add_task(
    add, name="add_calcfunction", x=wg.tasks.add_calcjob.outputs.sum, y=4
)  # We can access outputs as other tasks in WorkGraph
# (7*2)+3 = 17
wg.add_task(
    MultiplyAddWorkChain,
    name="multadd_workchain",
    x=wg.tasks.add_calcfunction.outputs.result,
    y=2,
    z=3,
    code=bash_code,
)
# 17+2+3 = 22
wg.add_task(
    add_more,
    name="add_more_workfunction",
    x=wg.tasks.multadd_workchain.outputs.result,
    y=2,
    z=3,
)
wg.run()
assert wg.tasks.add_more_workfunction.outputs.result.value == 22
print("Result:", wg.tasks.add_more_workfunction.outputs.result.value)

# %%
wg.to_html()

# %%
from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Calling WorkGraph within a WorkChain
# ====================================
# Although WorkGraph utilizes Calcfunction, Calcjob, and WorkChain to build workflows, it's not intended to replace the established WorkChain workflow system. As WorkGraph becomes more popular for creating workflows, there may arise scenarios where integrating a WorkGraph workflow within a WorkChain is necessary.
#
# To ensure system consistency, it's essential that WorkChains can seamlessly invoke WorkGraphs.
#
# This tutorial will demonstrate how to integrate a WorkGraph within a WorkChain.
#
# ## A Test WorkChain
#
# Let's construct a simple WorkChain that calls a WorkGraph. Consider the following key points:
#
# - **Dynamic Namespace**: Since we we will pass the whole WorkGraph data as input, the WorkChain must define a dynamic namespace to manage WorkGraph inputs.
# - **WorkGraph Submission**: Use the `WorkGraphEngine` class to submit a WorkGraph.
# - **Customization**: Users can specify a `call_link_label` to customize the link label for the process call.
# - **Output Handling**: It's important to understand the output namespace of the WorkGraph to retrieve results correctly.
#

# %%
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
        inputs = {'workgraph_data': self.inputs.workgraph,
                 'metadata': {'call_link_label': 'workgraph'}}
        process = self.submit(WorkGraphEngine, **inputs)
        self.to_context(workgraph=process)

    def results(self):
        # make sure that workgraph has outputs: `sum` and `product`
        self.out('sum', self.ctx.workgraph.outputs.sum)
        self.out('product', self.ctx.workgraph.outputs.product)
    


# %%
# Prepare the Inputs
# ------------------
# In this section, we'll create a basic `add_multiply` WorkGraph and pass it as input to the WorkChain.
#
# Key points to consider:
#
# - **Input Conversion**: Since WorkGraphs cannot be directly converted to AiiDA nodes, it is necessary to export the WorkGraph to a dictionary format before passing it to the WorkChain.
# - **Output Exposure**: Configure the WorkGraph to properly expose the outputs. This setup is crucial for the WorkChain to successfully retrieve the outputs. Note that the specific output names and the number of outputs will vary based on the particular use case.
#

# %%
from aiida.engine import calcfunction, run_get_node
from aiida import orm, load_profile
from aiida_workgraph import WorkGraph

load_profile()


@calcfunction
def add(x, y):
    return x + y

@calcfunction
def multiply(x, y):
    return x * y

# Create a `add_multiply` WorkGraph
wg = WorkGraph("add_multiply_workflow")
wg.add_task(add, name="add", x=orm.Int(1), y=orm.Int(2))
wg.add_task(multiply, name="multiply", x = wg.tasks.add.outputs.result,
                y=orm.Int(2))
# Define the outputs of the WorkGraph, which are exposed from the `multiply` and `add` tasks
#wg.outputs.product = wg.tasks.multiply.result
#wg.outputs.sum = wg.tasks.add.result
## We export the WorkGraph to a dictionary and pass it as input to the WorkChain
#inputs={"workgraph": wg.to_dict()}
#result, node = run_get_node(TestWorkChain, **inputs)

# %% [markdown]
# Print out the result and generate the node graph and visualize it.

# %%
from aiida_workgraph.utils import generate_node_graph
print("WorkChain results:")
#print("sum:    ", result['sum'])
#print("product:", result['product'])
#generate_node_graph(node.pk)

# %%
# Summary
# -------
# In this tutorial, we've covered how to integrate a dynamic WorkGraph within a WorkChain, addressing key aspects such as dynamic input management, process submission, and output retrieval.
#

# %%
# Further reading
# ===============
# One can also use ``WorkGraph`` inside a ``WorkChain``, please refer to the `Calling WorkGraph within a WorkChain <workchain_call_workgraph.ipynb>`_ for more details.
