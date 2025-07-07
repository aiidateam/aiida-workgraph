# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---
"""
=================================================
How to use aiida-core components inside WorkGraph
=================================================
"""

# %%
# Introduction
# ============
# If you’re already familiar with AiiDA, you may be interested in how to integrate ``aiida-core`` components—such as ``CalcJob``, ``calcfunction``, ``WorkChain``, and ``workfunction`` within a ``WorkGraph``. This integration enables you to seamlessly reuse existing workflows built with these paradigms inside a more flexible graph-based structure.
#
# Adding these components as tasks in an aiida-workgraph is straightforward, and their inputs and outputs can be easily connected to form complex workflows.
# .. note::
#
#    This guide assumes prior knowledge of AiiDA’s core components. If you’re unfamiliar with them, refer to the official documentation on `Calculations <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/index.html>`_ and `Workflows <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/workflows/index.html>`_.

# %%
# Preamble
# ========
# We need to initialize the computer and code resources to establish where the workflow is run.

from aiida import load_profile

load_profile()

from aiida import orm
from aiida.common.exceptions import NotExistent

try:
    computer = orm.Computer.collection.get(label="localhost")
except NotExistent:
    import tempfile

    computer = orm.Computer(
        label="localhost",
        description="localhoscomputer set up by test manager",
        hostname=label,
        workdir=str(tempfile.mkdtemp()),
        transport_type="core.local",
        scheduler_type="core.direct",
    )
    computer.store()
    computer.set_minimum_job_poll_interval(0.0)
    computer.configure()

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
# Integrating AiiDA's components
# ==============================
# We first have a look how the ``aiida-core`` components are represented in the workgraph widget as single tasks. We will see that the inputs and outputs that the components are visible in the widget.

from aiida_workgraph import WorkGraph

wg = WorkGraph("aiida_components")

# %%
# CalcJob
# -------
# As ``CalcJob`` we use the ``ArithmeticAddCalculation`` that adds to numbers ``x+y``.

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
# We create a ``workfunction`` that adds two three numbers ``x+y+z``.

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
# 1+2 = 3
wg.add_task(
    ArithmeticAddCalculation, name="add_calcjob", x=1, y=2, code=bash_code
)  # Note, that we add the code here
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

# %%
wg.to_html()

# %%
from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# Further reading
# ===============
# One can also use ``WorkGraph`` inside a ``WorkChain``, please refer to the [Calling WorkGraph within a WorkChain](workchain_call_workgraph.ipynb) for more details.
