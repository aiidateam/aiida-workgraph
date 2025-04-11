from aiida_workgraph.cli.cmd_workgraph import workgraph
import click
from aiida.cmdline.utils import decorators, echo
from aiida.cmdline.params import arguments, options
from aiida_workgraph.engine.scheduler.client import (
    get_scheduler_client,
    get_all_schedulers,
    get_scheduler,
)
from aiida import load_profile
from aiida_workgraph.engine.scheduler.scheduler import Scheduler
import kiwipy


@workgraph.group("scheduler")
def scheduler():
    """Commands to manage the scheduler process."""


@scheduler.command("list")
def scheduler_list():
    """Show a list of scheduler"""
    load_profile()
    schedulers = get_all_schedulers()
    echo.echo_formatted_list(schedulers, ["name"])


@scheduler.command("delete")
@click.argument("name", required=False, type=str)
@decorators.requires_loaded_profile()
def scheduler_delete(name):
    """Delete a scheduler."""
    from aiida.tools import delete_nodes

    scheduler = get_scheduler(name=name)
    status = Scheduler.get_status(name=scheduler.name)
    if status:
        echo.echo_error(
            f"Scheduler `{scheduler.name}` is running, please stop it first."
        )
        return

    _, was_deleted = delete_nodes([scheduler.pk], dry_run=False)
    if was_deleted:
        echo.echo_success(f"Scheduler `{name}` was deleted.")


@scheduler.command()
@click.argument("name", required=False, type=str)
@click.option("--foreground", is_flag=True, help="Run in foreground.")
@click.option(
    "--max-calcjob", type=int, required=False, help="Maximum number of calcjobs."
)
@click.option(
    "--max-process", type=int, required=False, help="Maximum number of processes."
)
@decorators.requires_broker
@decorators.check_circus_zmq_version
def start_scheduler(name, max_calcjob, max_process, foreground):
    """Start the scheduler application."""
    from aiida_workgraph.engine.scheduler.client import start_scheduler

    click.echo("Starting the scheduler process...")

    start_scheduler(
        name=name,
        max_calcjob=max_calcjob,
        max_process=max_process,
        foreground=foreground,
    )


@scheduler.command()
@click.argument("name", required=False, type=str)
@click.option("--foreground", is_flag=True, help="Run in foreground.")
@click.option(
    "--max-calcjob", type=int, required=False, help="Maximum number of calcjobs."
)
@click.option(
    "--max-process", type=int, required=False, help="Maximum number of processes."
)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@decorators.check_circus_zmq_version
def start(foreground, name, max_calcjob, max_process, timeout):
    """Start the scheduler application."""
    from aiida_workgraph.engine.scheduler.client import get_scheduler_client

    click.echo("Starting the scheduler ...")
    client = get_scheduler_client(scheduler_name=name)
    client.start_daemon(
        max_calcjob=max_calcjob,
        max_process=max_process,
        foreground=foreground,
        timeout=timeout,
    )


@scheduler.command()
@click.option("--no-wait", is_flag=True, help="Do not wait for confirmation.")
@click.argument("name", required=False, type=str)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@click.pass_context
def stop(ctx, name, no_wait, timeout):
    """Stop the scheduler daemon.

    Returns exit code 0 if the daemon was shut down successfully (or was not running), non-zero if there was an error.
    """
    if name:
        scheduler = get_scheduler(name=name)
        if not scheduler:
            echo.echo_error(f"Scheduler `{name}` not found.")
            ctx.exit(1)
        schedulers = [scheduler]
    else:
        schedulers = get_all_schedulers()

    for scheduler in schedulers:
        echo.echo("Scheduler: ", fg=echo.COLORS["report"], bold=True, nl=False)
        echo.echo(f"{scheduler.name}", bold=True)
        echo.echo("Stopping the scheduler... ")
        try:
            client = get_scheduler_client(scheduler_name=scheduler.name)
            client.stop_daemon(wait=not no_wait, timeout=timeout)
        except Exception as exception:
            echo.echo_error(f"Failed to stop the scheduler: {exception}")


@scheduler.command(hidden=True)
@click.argument(
    "name",
    required=False,
    type=str,
)
@click.option(
    "--max-calcjob", type=int, required=False, help="Maximum number of calcjobs."
)
@click.option(
    "--max-process", type=int, required=False, help="Maximum number of processes."
)
@click.option("--foreground", is_flag=True, help="Run in foreground.")
@decorators.with_dbenv()
@decorators.requires_broker
@decorators.check_circus_zmq_version
def start_circus(name, max_calcjob, max_process, foreground):
    """This will actually launch the circus daemon, either daemonized in the background or in the foreground.

    If run in the foreground all logs are redirected to stdout.

    .. note:: this should not be called directly from the commandline!
    """

    get_scheduler_client(scheduler_name=name)._start_daemon(
        max_calcjob=max_calcjob, max_process=max_process, foreground=foreground
    )


@scheduler.command()
@click.argument(
    "name",
    required=False,
    type=str,
)
@options.TIMEOUT(default=None, required=False, type=int)
@click.pass_context
@decorators.requires_loaded_profile()
@decorators.requires_broker
def status(ctx, name, timeout):
    """Print the status of the scheduler daemon.

    Returns exit code 0 if all requested daemons are running, else exit code 3.
    """
    from aiida.engine.daemon.client import DaemonException

    if name:
        schedulers = [get_scheduler(name=name)]
    else:
        schedulers = get_all_schedulers()

    print(
        "Name                status      pk   waiting  running   calcjob  max_calcjob max_process"
    )
    # for scheduler in schedulers:

    for scheduler in schedulers:
        echo.echo(f"{scheduler.name:<20}", bold=True, nl=False)

        try:
            result = Scheduler.get_status(name=scheduler.name)
            if result:
                echo.echo("Running", fg=echo.COLORS["success"], nl=False)
            else:
                echo.echo("Stopped", fg=echo.COLORS["error"], nl=False)
            echo.echo(f"  {scheduler.pk:<10}", nl=False)
            echo.echo(f"  {len(scheduler.waiting_process):<8}", nl=False)
            echo.echo(f"  {len(scheduler.running_process):<8}", nl=False)
            echo.echo(f"  {len(scheduler.running_calcjob):<8}", nl=False)
            echo.echo(f"  {scheduler.max_calcjob:<8}", nl=False)
            echo.echo(f"  {scheduler.max_process:<8}")
        except DaemonException as exception:
            echo.echo_error(str(exception))
            continue


@scheduler.command("show")
@click.argument(
    "name",
    required=False,
    type=str,
)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_loaded_profile()
@decorators.requires_broker
def scheduler_show(name, timeout):
    """Show details for a scheduler."""
    from aiida.tools.query.calculation import CalculationQueryBuilder
    from tabulate import tabulate

    builder = CalculationQueryBuilder()
    scheduler = get_scheduler(name=name)
    if not scheduler:
        echo.echo_error(f"Scheduler `{name}` not found.")
        return
    waiting_process = scheduler.waiting_process
    running_process = scheduler.running_process
    processes = waiting_process + running_process
    priorities = scheduler.get_process_priority()
    filters = None
    if processes:
        filters = {"id": {"in": processes}}

        query_set = builder.get_query_set(filters=filters)
        # project = default_projections()
        project = ("pk", "ctime", "process_label", "state")
        projected = builder.get_projected(query_set, projections=project)
        headers = projected.pop(0)
    else:
        projected = []
    # for all running processes
    for i, process in enumerate(projected):
        pk = process[0]
        # get the priority for the process
        if pk in priorities:
            process = list(process)
            process.append(priorities[pk])
            # process.append("Waiting")
            projected[i] = tuple(process)
        else:
            process = list(process)
            process.append(None)
            # process.append("Running")
            projected[i] = tuple(process)
    headers = ["PK", "Created", "Process label", "Process State", "Priorities"]

    echo.echo_report(f"Scheduler: {name}")
    tabulated = tabulate(projected, headers=headers)
    echo.echo(tabulated)
    echo.echo(f"\nTotal results: {len(projected)}\n")

    scheduler = get_scheduler(name=name)
    echo.echo(f"name: {scheduler.name}")
    echo.echo(f"pk: {scheduler.pk}")
    echo.echo(f"running_process: {len(running_process)}")
    echo.echo(f"waiting_process: {len(waiting_process)}")
    echo.echo(f"running_calcjob: {len(scheduler.running_calcjob)}")
    echo.echo(f"max_calcjob: {scheduler.max_calcjob}")
    echo.echo(f"max_process: {scheduler.max_process}")


@scheduler.command()
@click.argument("name", required=True, type=str)
@click.argument("max_calcjob", required=True, type=int)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@decorators.check_circus_zmq_version
def set_max_calcjob(name, max_calcjob, timeout):
    """Start the scheduler application."""

    try:
        Scheduler.set_max_calcjob(name=name, max_calcjob=max_calcjob)
    except kiwipy.exceptions.UnroutableError:
        echo.echo_error(
            f"Failed to set max_process for scheduler {name}. Is the scheduler running?"
        )


@scheduler.command()
@click.argument("name", required=True, type=str)
@click.argument("max_process", required=True, type=int)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@decorators.check_circus_zmq_version
def set_max_process(name, max_process, timeout):
    """Start the scheduler application."""
    try:
        Scheduler.set_max_process(name=name, max_process=max_process)
    except kiwipy.exceptions.UnroutableError:
        echo.echo_error(
            f"Failed to set max_process for scheduler {name}. Is the scheduler running?"
        )


@scheduler.command()
@click.argument("name", required=True, type=str)
@arguments.PROCESSES()
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@decorators.check_circus_zmq_version
def play_processes(name, processes, timeout):
    """Start the scheduler application."""
    pks = [p.pk for p in processes]
    try:
        Scheduler.play_processes(name=name, pks=pks, timeout=timeout)
    except kiwipy.exceptions.UnroutableError:
        echo.echo_error(
            f"Failed to set max_process for scheduler {name}. Is the scheduler running?"
        )
