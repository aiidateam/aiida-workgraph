from aiida.engine.daemon.client import DaemonClient
import shutil
from aiida.manage.manager import get_manager
from aiida.common.exceptions import ConfigurationError
import os

WORKGRAPH_BIN = shutil.which("workgraph")


class SchedulerClient(DaemonClient):
    """Client for interacting with the scheduler daemon."""

    _DAEMON_NAME = "scheduler-{name}"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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

        return {
            "circus": {
                "log": str(
                    DAEMON_LOG_DIR / f"circus-scheduler-{self.profile.name}.log"
                ),
                "pid": str(DAEMON_DIR / f"circus-scheduler-{self.profile.name}.pid"),
                "port": str(DAEMON_DIR / f"circus-scheduler-{self.profile.name}.port"),
                "socket": {
                    "file": str(
                        DAEMON_DIR / f"circus-scheduler-{self.profile.name}.sockets"
                    ),
                    "controller": "circus.c.sock",
                    "pubsub": "circus.p.sock",
                    "stats": "circus.s.sock",
                },
            },
            "daemon": {
                "log": str(DAEMON_LOG_DIR / f"aiida-{self.profile.name}.log"),
                "pid": str(DAEMON_DIR / f"aiida-{self.profile.name}.pid"),
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
        self, number_workers: int = 1, foreground: bool = False
    ) -> list[str]:
        """Return the command to start the daemon.

        :param number_workers: Number of daemon workers to start.
        :param foreground: Whether to launch the subprocess in the background or not.
        """
        command = [
            self._workgraph_bin,
            "-p",
            self.profile.name,
            "scheduler",
            "start-circus",
        ]

        if foreground:
            command.append("--foreground")

        return command

    @property
    def cmd_start_daemon_worker(self) -> list[str]:
        """Return the command to start a daemon worker process."""
        return [self._workgraph_bin, "-p", self.profile.name, "scheduler", "worker"]

    def _start_daemon(self, foreground: bool = False) -> None:
        """Start the daemon.

        .. warning:: This will daemonize the current process and put it in the background. It is most likely not what
            you want to call if you want to start the daemon from the Python API. Instead you probably will want to use
            the :meth:`aiida.engine.daemon.client.DaemonClient.start_daemon` function instead.

        :param number_workers: Number of daemon workers to start.
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
                    "cmd": " ".join(self.cmd_start_daemon_worker),
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


def get_scheduler_client(profile_name: str | None = None) -> "SchedulerClient":
    """Return the daemon client for the given profile or the current profile if not specified.

    :param profile_name: Optional profile name to load.
    :return: The daemon client.

    :raises aiida.common.MissingConfigurationError: if the configuration file cannot be found.
    :raises aiida.common.ProfileConfigurationError: if the given profile does not exist.
    """
    profile = get_manager().load_profile(profile_name)
    return SchedulerClient(profile)


def get_scheduler():
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
    qb.append(WorkGraphScheduler, filters=filters, project=projections, tag="process")
    results = qb.all()
    if len(results) == 0:
        raise ValueError("No scheduler found. Please start the scheduler first.")
    scheduler_id = results[0][0]
    return scheduler_id
