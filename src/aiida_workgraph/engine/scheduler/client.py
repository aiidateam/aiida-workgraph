from __future__ import annotations
from aiida.engine.daemon.client import (
    DaemonClient,
    DaemonException,
    DaemonTimeoutException,
)
import shutil
from aiida.manage.manager import get_manager
from aiida.common.exceptions import ConfigurationError
import os
from typing import Optional
from aiida.common.log import AIIDA_LOGGER
from typing import List
import subprocess
from aiida_workgraph.orm.scheduler import SchedulerNode

WORKGRAPH_BIN = shutil.which("workgraph")
LOGGER = AIIDA_LOGGER.getChild("engine.launch")


class SchedulerClient(DaemonClient):
    """Client for interacting with the scheduler daemon."""

    _DAEMON_NAME = "aiida-{profile_name}-{scheduler_name}"

    def __init__(self, scheduler_name, *args, **kwargs):
        self.scheduler_name = scheduler_name
        super().__init__(*args, **kwargs)

    @property
    def daemon_name(self) -> str:
        """Get the daemon name which is tied to the profile name."""
        return self._DAEMON_NAME.format(
            profile_name=self.profile.name, scheduler_name=self.scheduler_name
        )

    @property
    def _workgraph_bin(self) -> str:
        """Return the absolute path to the ``verdi`` binary.

        :raises ConfigurationError: If the path to ``verdi`` could not be found
        """
        if WORKGRAPH_BIN is None:
            raise ConfigurationError(
                "Unable to find 'verdi' in the path. Make sure that you are working "
                "in a virtual environment, or that at least the 'verdi' executable is on the PATH"
            )

        return WORKGRAPH_BIN

    @property
    def filepaths(self):
        """Return the filepaths used by this profile.

        :return: a dictionary of filepaths
        """
        from aiida.manage.configuration.settings import DAEMON_DIR, DAEMON_LOG_DIR

        # for aiida-core >= 2.7
        # from aiida.manage.configuration.settings import AiiDAConfigPathResolver
        # config = AiiDAConfigPathResolver()
        # DAEMON_LOG_DIR = config.daemon_log_dir
        # DAEMON_DIR = config.daemon_dir

        return {
            "circus": {
                "log": str(
                    DAEMON_LOG_DIR
                    / f"circus-scheduler-{self.profile.name}-{self.scheduler_name}.log"
                ),
                "pid": str(
                    DAEMON_DIR
                    / f"circus-scheduler-{self.profile.name}-{self.scheduler_name}.pid"
                ),
                "port": str(
                    DAEMON_DIR
                    / f"circus-scheduler-{self.profile.name}-{self.scheduler_name}.port"
                ),
                "socket": {
                    "file": str(
                        DAEMON_DIR
                        / f"circus-scheduler-{self.profile.name}-{self.scheduler_name}.sockets"
                    ),
                    "controller": "circus.c.sock",
                    "pubsub": "circus.p.sock",
                    "stats": "circus.s.sock",
                },
            },
            "daemon": {
                "log": str(
                    DAEMON_LOG_DIR
                    / f"aiida-scheduler-{self.profile.name}-{self.scheduler_name}.log"
                ),
                "pid": str(
                    DAEMON_DIR
                    / f"aiida-scheduler-{self.profile.name}-{self.scheduler_name}.pid"
                ),
            },
        }

    @property
    def circus_log_file(self) -> str:
        return self.filepaths["circus"]["log"]

    @property
    def circus_pid_file(self) -> str:
        return self.filepaths["circus"]["pid"]

    @property
    def circus_port_file(self) -> str:
        return self.filepaths["circus"]["port"]

    @property
    def circus_socket_file(self) -> str:
        return self.filepaths["circus"]["socket"]["file"]

    @property
    def circus_socket_endpoints(self) -> dict[str, str]:
        return self.filepaths["circus"]["socket"]

    @property
    def daemon_log_file(self) -> str:
        return self.filepaths["daemon"]["log"]

    @property
    def daemon_pid_file(self) -> str:
        return self.filepaths["daemon"]["pid"]

    def cmd_start_daemon(
        self,
        max_calcjobs: int | None = None,
        max_workflows: int | None = None,
        max_processes: int | None = None,
        foreground: bool = False,
    ) -> list[str]:
        """Return the command to start the daemon.

        :param foreground: Whether to launch the subprocess in the background or not.
        """
        command = [
            self._workgraph_bin,
            "-p",
            self.profile.name,
            "scheduler",
            "start-circus",
            self.scheduler_name,
        ]

        if max_calcjobs is not None:
            command.append(f"--max-calcjobs={max_calcjobs}")
        if max_workflows is not None:
            command.append(f"--max-workflows={max_workflows}")
        if max_processes is not None:
            command.append(f"--max-processes={max_processes}")
        if foreground:
            command.append("--foreground")

        return command

    def cmd_start_daemon_worker(
        self,
        max_calcjobs: int | None = None,
        max_workflows: int | None = None,
        max_processes: int | None = None,
    ) -> list[str]:
        """Return the command to start a daemon worker process."""
        command = [
            self._workgraph_bin,
            "-p",
            self.profile.name,
            "scheduler",
            "start-scheduler",
            self.scheduler_name,
        ]

        if max_calcjobs is not None:
            command.append(f"--max-calcjobs={max_calcjobs}")
        if max_workflows is not None:
            command.append(f"--max-workflows={max_workflows}")
        if max_processes is not None:
            command.append(f"--max-processes={max_processes}")

        return command

    def get_scheduler_node(self) -> Optional[SchedulerNode]:
        """Return the scheduler node with the given name.

        :return: The scheduler node or None if not found.
        """
        from aiida import orm

        qb = orm.QueryBuilder()
        qb.append(SchedulerNode, filters={"attributes.name": self.scheduler_name})
        return qb.first()[0] if qb.count() else None

    def start_daemon(
        self,
        max_calcjobs: int | None = None,
        max_workflows: int | None = None,
        max_processes: int | None = None,
        foreground: bool = False,
        wait: bool = True,
        timeout: int | None = None,
    ) -> None:
        """Start the daemon in a sub process running in the background.

        :param name: Name of the scheduler to start.
        :param foreground: Whether to launch the subprocess in the background or not.
        :param wait: Boolean to indicate whether to wait for the result of the command.
        :param timeout: Optional timeout to set for trying to reach the circus daemon after the subprocess has started.
            Default is set on the client upon instantiation taken from the ``daemon.timeout`` config option.
        :raises DaemonException: If the command to start the daemon subprocess excepts.
        :raises DaemonTimeoutException: If the daemon starts but then is unresponsive or in an unexpected state.
        """
        self._clean_potentially_stale_pid_file()

        env = self.get_env()
        command = self.cmd_start_daemon(
            max_calcjobs=max_calcjobs,
            max_workflows=max_workflows,
            max_processes=max_processes,
            foreground=foreground,
        )
        timeout = timeout or self._daemon_timeout

        try:
            subprocess.check_output(command, env=env, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exception:
            raise DaemonException("The daemon failed to start.") from exception

        if not wait:
            return

        self._await_condition(
            lambda: self.is_daemon_running,
            DaemonTimeoutException(
                f"The daemon failed to start or is unresponsive after {timeout} seconds."
            ),
            timeout=timeout,
        )

    def _start_daemon(
        self,
        max_calcjobs: int | None = None,
        max_workflows: int | None = None,
        max_processes: int | None = None,
        foreground: bool = False,
    ) -> None:
        """Start the daemon.

        .. warning:: This will daemonize the current process and put it in the background. It is most likely not what
            you want to call if you want to start the daemon from the Python API. Instead you probably will want to use
            the :meth:`aiida.engine.daemon.client.DaemonClient.start_daemon` function instead.

        :param name: Name of the scheduler to start.
        :param foreground: Whether to launch the subprocess in the background or not.
        """
        from circus import get_arbiter
        from circus import logger as circus_logger
        from circus.circusd import daemonize
        from circus.pidfile import Pidfile
        from circus.util import check_future_exception_and_log, configure_logger

        loglevel = self.loglevel
        logoutput = "-"

        if not foreground:
            logoutput = self.circus_log_file

        arbiter_config = {
            "controller": self.get_controller_endpoint(),
            "pubsub_endpoint": self.get_pubsub_endpoint(),
            "stats_endpoint": self.get_stats_endpoint(),
            "logoutput": logoutput,
            "loglevel": loglevel,
            "debug": False,
            "statsd": True,
            "pidfile": self.circus_pid_file,
            "watchers": [
                {
                    "cmd": " ".join(
                        self.cmd_start_daemon_worker(
                            max_calcjobs=max_calcjobs,
                            max_workflows=max_workflows,
                            max_processes=max_processes,
                        )
                    ),
                    "name": self.daemon_name,
                    "numprocesses": 1,
                    "virtualenv": self.virtualenv,
                    "copy_env": True,
                    "stdout_stream": {
                        "class": "FileStream",
                        "filename": self.daemon_log_file,
                    },
                    "stderr_stream": {
                        "class": "FileStream",
                        "filename": self.daemon_log_file,
                    },
                    "env": self.get_env(),
                }
            ],
        }

        if not foreground:
            daemonize()

        arbiter = get_arbiter(**arbiter_config)
        pidfile = Pidfile(arbiter.pidfile)
        pidfile.create(os.getpid())

        # Configure the logger
        loggerconfig = arbiter.loggerconfig or None
        configure_logger(circus_logger, loglevel, logoutput, loggerconfig)

        # Main loop
        should_restart = True

        while should_restart:
            try:
                future = arbiter.start()
                should_restart = False
                if check_future_exception_and_log(future) is None:
                    should_restart = arbiter._restarting
            except Exception as exception:
                # Emergency stop
                arbiter.loop.run_sync(arbiter._emergency_stop)
                raise exception
            except KeyboardInterrupt:
                pass
            finally:
                arbiter = None
                if pidfile is not None:
                    pidfile.unlink()


def get_scheduler_client(
    scheduler_name: str, profile_name: Optional[str] = None
) -> "SchedulerClient":
    """Return the daemon client for the given profile or the current profile if not specified.

    :param profile_name: Optional profile name to load.
    :return: The daemon client.

    :raises aiida.common.MissingConfigurationError: if the configuration file cannot be found.
    :raises aiida.common.ProfileConfigurationError: if the given profile does not exist.
    """
    profile = get_manager().load_profile(profile_name)
    return SchedulerClient(scheduler_name=scheduler_name, profile=profile)


def get_all_scheduler_nodes() -> List[str]:
    from aiida import orm

    qb = orm.QueryBuilder()
    qb.append(SchedulerNode)
    return [scheduler[0] for scheduler in qb.all()]


def get_scheduler_node(name: str) -> Optional[SchedulerNode]:
    """Return the scheduler node with the given name.

    :param name: The name of the scheduler.
    :return: The scheduler node or None if not found.
    """
    from aiida import orm

    qb = orm.QueryBuilder()
    qb.append(SchedulerNode, filters={"attributes.name": name})
    return qb.first()[0] if qb.count() else None


def start_scheduler(
    name: str,
    max_calcjobs: int | None = None,
    max_workflows: int | None = None,
    max_processes: int | None = None,
    foreground: bool = False,
) -> None:
    """Start a scheduler worker for the currently configured profile.

    :param foreground: If true, the logging will be configured to write to stdout, otherwise it will be configured to
        write to the scheduler log file.
    """
    import asyncio
    import signal
    import sys
    from aiida_workgraph.engine.scheduler.client import get_scheduler_client
    from aiida_workgraph.engine.scheduler.scheduler import Scheduler

    from aiida.common.log import configure_logging
    from aiida.manage import get_config_option
    from aiida.engine.daemon.worker import shutdown_worker

    daemon_client = get_scheduler_client(scheduler_name=name)
    configure_logging(
        daemon=not foreground, daemon_log_file=daemon_client.daemon_log_file
    )

    LOGGER.debug(f"sys.executable: {sys.executable}")
    LOGGER.debug(f"sys.path: {sys.path}")

    try:
        scheduler = Scheduler(
            name=name,
            max_calcjobs=max_calcjobs,
            max_workflows=max_workflows,
            max_processes=max_processes,
        )
    except Exception:
        LOGGER.exception("daemon worker failed to start")
        raise

    if isinstance(rlimit := get_config_option("daemon.recursion_limit"), int):
        LOGGER.info("Setting maximum recursion limit of daemon worker to %s", rlimit)
        sys.setrecursionlimit(rlimit)

    signals = (signal.SIGTERM, signal.SIGINT)
    for s in signals:
        # https://github.com/python/mypy/issues/12557
        scheduler._loop.add_signal_handler(
            s, lambda s=s: asyncio.create_task(shutdown_worker(scheduler))
        )

    try:
        LOGGER.info(f"Starting scheduler: {name}")
        scheduler.start()
    except SystemError as exception:
        LOGGER.info("Received a SystemError: %s", exception)
        scheduler.close()

    LOGGER.info(f"Scheduler: {name} started")
