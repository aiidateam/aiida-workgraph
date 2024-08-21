"""`verdi process` command."""
import click

from aiida_workgraph.cli.cmd_workgraph import workgraph
from aiida.cmdline.params import arguments, options
from aiida.cmdline.utils import decorators, echo
from aiida_workgraph.cli.query_workgraph import WorkGraphQueryBuilder

REPAIR_INSTRUCTIONS = """\
If one ore more processes are unreachable, you can run the following commands to try and repair them:

    verdi daemon stop
    verdi process repair
    verdi daemon start
"""


def default_projections():
    """Return list of default projections for the ``--project`` option of ``verdi process list``.

    This indirection is necessary to prevent loading the imported module which slows down tab-completion.
    """
    return WorkGraphQueryBuilder.default_projections


@workgraph.group("task")
def workgraph_task():
    """Inspect and manage processes."""


@workgraph_task.command("pause")
@arguments.PROCESS()
@click.argument("tasks", nargs=-1)
@options.TIMEOUT()
@options.WAIT()
@decorators.with_dbenv()
def task_pause(process, tasks, timeout, wait):
    """Pause task."""
    from aiida.engine.processes import control
    from aiida_workgraph.utils.control import pause_tasks

    try:
        _, msg = pause_tasks(process.pk, tasks, timeout, wait)
    except control.ProcessTimeoutException as exception:
        echo.echo_critical(f"{exception}\n{REPAIR_INSTRUCTIONS}")


@workgraph_task.command("play")
@arguments.PROCESS()
@click.argument("tasks", nargs=-1)
@options.TIMEOUT()
@options.WAIT()
@decorators.with_dbenv()
def task_play(process, tasks, timeout, wait):
    """Play task."""
    from aiida.engine.processes import control
    from aiida_workgraph.utils.control import play_tasks

    try:
        _, msg = play_tasks(process.pk, tasks, timeout, wait)
    except control.ProcessTimeoutException as exception:
        echo.echo_critical(f"{exception}\n{REPAIR_INSTRUCTIONS}")


@workgraph_task.command("skip")
@arguments.PROCESS()
@click.argument("tasks", nargs=-1)
@options.TIMEOUT()
@options.WAIT()
@decorators.with_dbenv()
def task_skip(process, tasks, timeout, wait):
    """Skip task."""
    from aiida.engine.processes import control
    from aiida_workgraph.utils.control import skip_tasks

    for task in tasks:
        try:
            skip_tasks(process.pk, task, timeout, wait)

        except control.ProcessTimeoutException as exception:
            echo.echo_critical(f"{exception}\n{REPAIR_INSTRUCTIONS}")


@workgraph_task.command("kill")
@arguments.PROCESS()
@click.argument("tasks", nargs=-1)
@options.TIMEOUT()
@options.WAIT()
@decorators.with_dbenv()
def task_kill(process, tasks, timeout, wait):
    """Kill task."""
    from aiida.engine.processes import control
    from aiida_workgraph.utils.control import kill_tasks

    print("tasks", tasks)
    try:
        kill_tasks(process.pk, tasks, timeout, wait)
    except control.ProcessTimeoutException as exception:
        echo.echo_critical(f"{exception}\n{REPAIR_INSTRUCTIONS}")
