from aiida_workgraph.cli.cmd_workgraph import workgraph
import click
from aiida.cmdline.utils import decorators, echo
from aiida.cmdline.commands.cmd_daemon import validate_daemon_workers
from aiida.cmdline.params import options
from aiida_workgraph.engine.scheduler.client import get_scheduler_client
import sys


@workgraph.group("scheduler")
def scheduler():
    """Commands to manage the scheduler process."""


@scheduler.command()
def worker():
    """Start the scheduler application."""
    from aiida_workgraph.engine.scheduler.client import start_scheduler_worker

    click.echo("Starting the scheduler worker...")

    start_scheduler_worker()


@scheduler.command()
@click.option("--foreground", is_flag=True, help="Run in foreground.")
@click.argument("number", required=False, type=int, callback=validate_daemon_workers)
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.with_dbenv()
@decorators.requires_broker
@decorators.check_circus_zmq_version
def start(foreground, number, timeout):
    """Start the scheduler application."""
    from aiida_workgraph.engine.scheduler.client import start_scheduler_process

    click.echo("Starting the scheduler process...")

    client = get_scheduler_client()
    client.start_daemon(number_workers=number, foreground=foreground, timeout=timeout)
    start_scheduler_process(number)


@scheduler.command()
@click.option("--no-wait", is_flag=True, help="Do not wait for confirmation.")
@click.option("--all", "all_profiles", is_flag=True, help="Stop all daemons.")
@options.TIMEOUT(default=None, required=False, type=int)
@decorators.requires_broker
@click.pass_context
def stop(ctx, no_wait, all_profiles, timeout):
    """Stop the scheduler daemon.

    Returns exit code 0 if the daemon was shut down successfully (or was not running), non-zero if there was an error.
    """
    if all_profiles is True:
        profiles = [
            profile
            for profile in ctx.obj.config.profiles
            if not profile.is_test_profile
        ]
    else:
        profiles = [ctx.obj.profile]

    for profile in profiles:
        echo.echo("Profile: ", fg=echo.COLORS["report"], bold=True, nl=False)
        echo.echo(f"{profile.name}", bold=True)
        echo.echo("Stopping the daemon... ", nl=False)
        try:
            client = get_scheduler_client()
            client.stop_daemon(wait=not no_wait, timeout=timeout)
        except Exception as exception:
            echo.echo_error(f"Failed to stop the daemon: {exception}")


@scheduler.command(hidden=True)
@click.option("--foreground", is_flag=True, help="Run in foreground.")
@click.argument("number", required=False, type=int, callback=validate_daemon_workers)
@decorators.with_dbenv()
@decorators.requires_broker
@decorators.check_circus_zmq_version
def start_circus(foreground, number):
    """This will actually launch the circus daemon, either daemonized in the background or in the foreground.

    If run in the foreground all logs are redirected to stdout.

    .. note:: this should not be called directly from the commandline!
    """

    get_scheduler_client()._start_daemon(number_workers=number, foreground=foreground)


@scheduler.command()
@click.option("--all", "all_profiles", is_flag=True, help="Show status of all daemons.")
@options.TIMEOUT(default=None, required=False, type=int)
@click.pass_context
@decorators.requires_loaded_profile()
@decorators.requires_broker
def status(ctx, all_profiles, timeout):
    """Print the status of the scheduler daemon.

    Returns exit code 0 if all requested daemons are running, else exit code 3.
    """
    from tabulate import tabulate

    from aiida.cmdline.utils.common import format_local_time
    from aiida.engine.daemon.client import DaemonException

    if all_profiles is True:
        profiles = [
            profile
            for profile in ctx.obj.config.profiles
            if not profile.is_test_profile
        ]
    else:
        profiles = [ctx.obj.profile]

    daemons_running = []

    for profile in profiles:
        client = get_scheduler_client(profile.name)
        echo.echo("Profile: ", fg=echo.COLORS["report"], bold=True, nl=False)
        echo.echo(f"{profile.name}", bold=True)

        try:
            client.get_status(timeout=timeout)
        except DaemonException as exception:
            echo.echo_error(str(exception))
            daemons_running.append(False)
            continue

        worker_response = client.get_worker_info()
        daemon_response = client.get_daemon_info()

        workers = []
        for pid, info in worker_response["info"].items():
            if isinstance(info, dict):
                row = [
                    pid,
                    info["mem"],
                    info["cpu"],
                    format_local_time(info["create_time"]),
                ]
            else:
                row = [pid, "-", "-", "-"]
            workers.append(row)

        if workers:
            workers_info = tabulate(
                workers, headers=["PID", "MEM %", "CPU %", "started"], tablefmt="simple"
            )
        else:
            workers_info = (
                "--> No workers are running. Use `verdi daemon incr` to start some!\n"
            )

        start_time = format_local_time(daemon_response["info"]["create_time"])
        echo.echo(
            f'Daemon is running as PID {daemon_response["info"]["pid"]} since {start_time}\n'
            f"Active workers [{len(workers)}]:\n{workers_info}\n"
            "Use `verdi daemon [incr | decr] [num]` to increase / decrease the number of workers"
        )

    if not all(daemons_running):
        sys.exit(3)
