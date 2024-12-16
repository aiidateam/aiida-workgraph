"""
=======================
ShellJob
=======================

"""

# %%
# Introduction
# ============
# ``ShellJob`` is a built-in task, which uses the `aiida-shell <https://aiida-shell.readthedocs.io/en/latest/>`_ package to run shell commands easily. Run any shell executable without writing a dedicated plugin or parser.
#
# This tutorial is based on the `docs <https://aiida-shell.readthedocs.io/en/latest/howto.html#>`_ of the ``aiida-shell``.
#
# Load the AiiDA profile.

from aiida import load_profile

load_profile()

# %%
# Running a shell command
# ========================
# Run a shell command without any arguments. Here we run the `date` command to show the date.
#

from aiida_workgraph import WorkGraph

wg = WorkGraph(name="test_shell_date")
date_task = wg.add_task("ShellJob", command="date")
wg.submit(wait=True)

# Print out the result:
print("\nResult: ", date_task.outputs.stdout.value.get_content())

# %%
# Under the hood, an AiiDA ``Code`` instance ``date`` will be created on the ``localhost`` computer. In addition, it is also
# possible to specify a different, remote computer. For the sake of this example, we create a mock remote computer (you
# would usually have this pre-configured already, e.g., your favorite HPC resource):

from aiida import orm

created, mock_remote_computer = orm.Computer.collection.get_or_create(
    label="my-remote-computer",
    description="A mock remote computer",
    hostname="my-remote-computer",
    workdir="/tmp/aiida",
    transport_type="core.ssh",
    scheduler_type="core.direct",
)
if created:
    mock_remote_computer.store()

# We can then specify the remote computer via the ``metadata``:

wg = WorkGraph(name="test_shell_date_remote")
date_task_remote = wg.add_task(
    "ShellJob",
    command="date",
    metadata={"computer": orm.load_computer("my-remote-computer")},
)

#%%
# In addition, it is possible to pass a pre-configured ``Code``. Let's again create a mock ``Code``:

from aiida_workgraph.utils import get_or_create_code

mock_remote_code = get_or_create_code(
    code_label="remote-date",
    computer="my-remote-computer",
    code_path="/usr/bin/date",
)

# Which we can now directly pass to the ``command`` argument of ``add_task``:

wg = WorkGraph(name="test_shell_date_preconfigured")
preconf_task = wg.add_task(
    "ShellJob", command=orm.load_code("remote-date@my-remote-computer")
)

# Note, we are not running or submitting the ``WorkGraph`` in the above two examples, as we are using mocked
# ``Computer`` and ``Code`` entities for demonstration purposes.

# %%
# Running a shell command with arguments
# =======================================
# To pass arguments to the shell command, pass them as a list of strings to the arguments keyword:
#


wg = WorkGraph(name="test_shell_date_with_arguments")
date_task = wg.add_task("ShellJob", command="date", arguments=["--iso-8601"])
wg.submit(wait=True)

# Print out the result:
print("\nResult: ", date_task.outputs.stdout.value.get_content())

# %%
# Running a shell command with files as arguments
# ===============================================
# For commands that take arguments that refer to files, pass those files using the nodes keyword. The keyword takes a dictionary of SinglefileData nodes. To specify where on the command line the files should be passed, use placeholder strings in the arguments keyword.
#


from aiida.orm import SinglefileData

wg = WorkGraph(name="test_shell_cat_with_file_arguments")
cat_task = wg.add_task(
    "ShellJob",
    command="cat",
    arguments=["{file_a}", "{file_b}"],
    nodes={
        "file_a": SinglefileData.from_string("string a"),
        "file_b": SinglefileData.from_string("string b"),
    },
)
wg.submit(wait=True)

# Print out the result:
print("\nResult: ", cat_task.outputs.stdout.value.get_content())

# %%
# Create a workflow
# =================
# We want to calculate `(x+y)*z` in two steps using the `expr` command. Each step will invole one `ShellJob`.
#
# If one wanted to run this workflow in AiiDA using CalcJob and WorkChain, one would have to write plugins for `expr` command, and a WorkChain to handle the workflow. With aiida-workgraph, this can be run with the following workgraph:
#

from aiida_workgraph import WorkGraph
from aiida.orm import Int
from aiida_shell.data import PickledData
from aiida import load_profile

load_profile()


def parser(dirpath):
    from aiida.orm import Int

    return {"result": Int((dirpath / "stdout").read_text().strip())}


# Create a workgraph
wg = WorkGraph(name="shell_add_mutiply_workflow")
expr_1 = wg.add_task(
    "ShellJob",
    name="expr_1",
    command="expr",
    arguments=["{x}", "+", "{y}"],
    nodes={"x": Int(2), "y": Int(3)},
    parser=PickledData(parser),
    parser_outputs=[{"name": "result"}],
)
expr_2 = wg.add_task(
    "ShellJob",
    name="expr_2",
    command="expr",
    arguments=["{result}", "*", "{z}"],
    nodes={"z": Int(4), "result": expr_1.outputs.result},
    parser=PickledData(parser),
    parser_outputs=[{"name": "result"}],
)
wg.to_html()


# %%
# Submit the workgraph
#

wg.submit(wait=True, timeout=100)
print("State of WorkGraph    : {}".format(wg.state))
print("Result               : {}".format(expr_2.node.outputs.result.value))


# %%
# Generate node graph from the AiiDA process
#

from aiida_workgraph.utils import generate_node_graph

generate_node_graph(wg.pk)

# %%
# What's Next
# ===========
# Overall, WorkGraph's ``ShellJob`` builds on top of ``aiida-shell`` and mirrors its functionality. So features
# implemented in ``aiida-shell`` should also be available for the ``ShellJob``. For more examples of ``aiida-shell``, please refer to its `docs <https://aiida-shell.readthedocs.io/en/latest/howto.html#>`_.
#
