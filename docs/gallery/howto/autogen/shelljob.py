"""
==============================
Run shell commands as a task
==============================

"""

# %%
# Introduction
# ------------
#
# This document demonstrates how to use the built-in ``shelljob`` task within AiiDA-WorkGraph, which leverages the `aiida-shell <https://aiida-shell.readthedocs.io/en/latest/>`_ package, to easily execute shell commands.
# This eliminates the need to write custom plugins or parsers for simple shell operations.
#
# This tutorial is based on the `docs <https://aiida-shell.readthedocs.io/en/latest/howto.html#>`_ of the `aiida-shell` package.

import typing as t

from aiida import load_profile, orm

from aiida_workgraph import dynamic, shelljob, task
from aiida_workgraph.utils import get_or_create_code

load_profile()

# %%
# Run a shell command
# -------------------
#
# To execute a shell command, you can use ``shelljob``.
# Here's an example that runs the `date` command to retrieve the current date:


@task.graph
def ShellDate() -> t.Annotated[dict, dynamic(t.Any)]:
    return shelljob(command='date')


wg = ShellDate.build()
wg.run()

print('\nResult: ', wg.outputs.stdout.value.get_content())

# %%
# .. note::
#
#    In the current version, it is required to annotate the return type of a task which returns the outputs of a ``shelljob`` with ``t.Annotated[dict, dynamic(t.Any)]``.
#    The annotation overrides the default ``result`` socket with that of the ``shelljob``.
#    In a future release, this will be handled via explicit exposing of the ``shelljob`` outputs.
#
# Under the hood, an AiiDA ``Code`` instance referencing the ``date`` shell command will be created on the ``localhost`` computer.
# In addition, it is also possible to specify a different remote computer.
# For the sake of this example, we create a mock remote computer (you would usually have this pre-configured already, e.g., your favorite HPC resource).


created, mock_remote_computer = orm.Computer.collection.get_or_create(
    label='my-remote-computer',
    description='A mock remote computer',
    hostname='my-remote-computer',
    workdir='/tmp/aiida',
    transport_type='core.ssh',
    scheduler_type='core.direct',
)
if created:
    mock_remote_computer.store()

# %%
# We can then specify the remote computer via the ``metadata``.


@task.graph
def RemoteShellDate(computer: str) -> t.Annotated[dict, dynamic(t.Any)]:
    return shelljob(
        command='date',
        metadata={'computer': orm.load_computer(computer)},
    )


# %%
# In addition, it is possible to pass a pre-configured ``Code``.
# Let's again create a mock ``Code``.


mock_remote_code = get_or_create_code(
    code_label='remote-date',
    computer='my-remote-computer',
    code_path='/usr/bin/date',
)

# %%
# We can now directly pass to the ``command`` argument of ``add_task``.


@task.graph
def RemoteShellDateWithCode(code: str) -> t.Annotated[dict, dynamic(t.Any)]:
    return shelljob(command=orm.load_code(code))


# %%
# Command arguments
# -----------------
#
# To pass arguments to the shell command, pass them as a list of strings to the arguments keyword.


@task.graph
def ShellDateWithArguments() -> t.Annotated[dict, dynamic(t.Any)]:
    return shelljob(command='date', arguments=['--iso-8601'])


wg = ShellDateWithArguments.build()
wg.run()

print('\nResult: ', wg.outputs.stdout.value.get_content())

# %%
# File arguments
# --------------
#
# For commands that take arguments referring to files, you pass those files using the nodes keyword.
# The keyword takes a dictionary of `SinglefileData` nodes.
# To specify where on the command line the files should be passed, use placeholder strings in the arguments keyword.


@task.graph
def ShellCatWithFileArguments() -> t.Annotated[dict, dynamic(t.Any)]:
    return shelljob(
        command='cat',
        arguments=['{file_a}', '{file_b}'],
        nodes={
            'file_a': orm.SinglefileData.from_string('string a'),
            'file_b': orm.SinglefileData.from_string('string b'),
        },
    )


wg = ShellCatWithFileArguments.build()
wg.run()

print('\nResult: ', wg.outputs.stdout.value.get_content())

# %%
# An example workflow: ``(x + y) * z``
# ------------------------------------
#
# We want to calculate ``(x + y) * z`` in two steps using the ``expr`` command.
# Each step will involve one `ShellJob`.
# To facilitate access to the result of each step, we define a custom parser that reads the output from the `stdout` file.


def parser(dirpath):
    return {'result': orm.Int((dirpath / 'stdout').read_text().strip())}


@task.graph
def ShellAddMultiply(x: int, y: int, z: int) -> int:
    the_sum = shelljob(
        command='expr',
        arguments=['{x}', '+', '{y}'],
        nodes={'x': x, 'y': y},
        parser=parser,
        parser_outputs=['result'],
    ).result
    the_product = shelljob(
        command='expr',
        arguments=['{x}', '*', '{y}'],
        nodes={'x': the_sum, 'y': z},
        parser=parser,
        parser_outputs=['result'],
    ).result
    return the_product


# Create a workgraph
wg = ShellAddMultiply.build(2, 3, 4)
wg.to_html()

# %%

wg.run()

print('State` of WorkGraph    : {}'.format(wg.state))
print('Result                : {}'.format(wg.outputs.result.value))
assert wg.outputs.result.value == 20


# %%
# Let's have a look at the provenance graph.

wg.generate_provenance_graph()

# %%
# What's Next
# -----------
#
# Overall, WorkGraph's `ShellJob` builds on top of ``aiida-shell`` and mirrors its functionality.
# Features implemented in ``aiida-shell`` are also be available for the ``ShellJob``.
# For more examples of ``aiida-shell``, please refer to its `docs <https://aiida-shell.readthedocs.io/en/latest/howto.html#>`_.
