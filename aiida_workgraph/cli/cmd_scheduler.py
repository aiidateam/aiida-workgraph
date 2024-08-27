from aiida_workgraph.cli.cmd_workgraph import workgraph
from aiida import orm
import click
import os
from pathlib import Path
from aiida.cmdline.utils import echo
from .cmd_graph import REPAIR_INSTRUCTIONS


REACT_PORT = "3000"


def get_package_root():
    """Returns the root directory of the package."""
    current_file = Path(__file__)
    # Root directory of your package
    return current_file.parent


def get_pid_file_path():
    """Get the path to the PID file in the desired directory."""
    from aiida.manage.configuration.settings import AIIDA_CONFIG_FOLDER

    return AIIDA_CONFIG_FOLDER / "scheduler_processes.pid"


@workgraph.group("scheduler")
def scheduler():
    """Commands to manage the scheduler process."""


@scheduler.command()
def start_worker():
    """Start the scheduler application."""
    from aiida_workgraph.engine.launch import start_scheduler_worker

    click.echo("Starting the scheduler worker...")

    start_scheduler_worker()


@scheduler.command()
def start():
    """Start the scheduler application."""
    from aiida_workgraph.engine.scheduler import WorkGraphScheduler
    from aiida.engine import submit

    click.echo("Starting the scheduler process...")

    pid_file_path = get_pid_file_path()
    # if the PID file already exists, check if the process is running
    if pid_file_path.exists():
        with open(pid_file_path, "r") as pid_file:
            for line in pid_file:
                _, pid = line.strip().split(":")
                if pid:
                    try:
                        node = orm.load_node(pid)
                        if node.is_sealed:
                            click.echo(
                                "PID file exists but no running scheduler process found."
                            )
                        else:
                            click.echo(
                                f"Scheduler process with PID {node.pk} already running."
                            )
                            return
                    except Exception:
                        click.echo(
                            "PID file exists but no running scheduler process found."
                        )

    with open(pid_file_path, "w") as pid_file:
        node = submit(WorkGraphScheduler)
        pid_file.write(f"Scheduler:{node.pk}\n")
        click.echo(f"Scheduler process started with PID {node.pk}.")


@scheduler.command()
def stop():
    """Stop the scheduler application."""
    from aiida.engine.processes import control

    pid_file_path = get_pid_file_path()

    if not pid_file_path.exists():
        click.echo("No running scheduler application found.")
        return

    with open(pid_file_path, "r") as pid_file:
        for line in pid_file:
            _, pid = line.strip().split(":")
            if pid:
                click.confirm(
                    "Are you sure you want to kill the scheduler process?", abort=True
                )
            process = orm.load_node(pid)
            try:
                message = "Killed through `verdi process kill`"
                control.kill_processes(
                    [process],
                    timeout=5,
                    wait=True,
                    message=message,
                )
            except control.ProcessTimeoutException as exception:
                echo.echo_critical(f"{exception}\n{REPAIR_INSTRUCTIONS}")
    os.remove(pid_file_path)


@scheduler.command()
def status():
    """Check the status of the scheduler application."""
    from aiida.orm import QueryBuilder
    from aiida_workgraph.engine.scheduler import WorkGraphScheduler

    qb = QueryBuilder()
    projections = ["id"]
    filters = {
        "or": [
            {"attributes.sealed": False},
            {"attributes": {"!has_key": "sealed"}},
        ]
    }
    qb.append(
        WorkGraphScheduler,
        filters=filters,
        project=projections,
        tag="process",
    )
    results = qb.all()
    if len(results) == 0:
        click.echo("No scheduler found. Please start the scheduler first.")
    else:
        click.echo(f"Scheduler process is running with PID: {results[0][0]}")
