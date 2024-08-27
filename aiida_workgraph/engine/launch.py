from __future__ import annotations

import time
import typing as t

from aiida.common import InvalidOperation
from aiida.common.log import AIIDA_LOGGER
from aiida.manage import manager
from aiida.orm import ProcessNode

from aiida.engine.processes.builder import ProcessBuilder
from aiida.engine.processes.functions import get_stack_size
from aiida.engine.processes.process import Process
from aiida.engine.utils import prepare_inputs, is_process_function

import signal
import sys
import inspect
from typing import (
    Type,
    Union,
)


from aiida.manage import get_manager

__all__ = ("run_get_node", "submit")

TYPE_RUN_PROCESS = t.Union[Process, t.Type[Process], ProcessBuilder]
# run can also be process function, but it is not clear what type this should be
TYPE_SUBMIT_PROCESS = t.Union[Process, t.Type[Process], ProcessBuilder]
LOGGER = AIIDA_LOGGER.getChild("engine.launch")


def run_get_node(
    process_class, *args, **kwargs
) -> tuple[dict[str, t.Any] | None, "ProcessNode"]:
    """Run the FunctionProcess with the supplied inputs in a local runner.
    :param args: input arguments to construct the FunctionProcess
    :param kwargs: input keyword arguments to construct the FunctionProcess
    :return: tuple of the outputs of the process and the process node
    """
    parent_pid = kwargs.pop("_parent_pid", None)
    frame_delta = 1000
    frame_count = get_stack_size()
    stack_limit = sys.getrecursionlimit()
    LOGGER.info(
        "Executing process function, current stack status: %d frames of %d",
        frame_count,
        stack_limit,
    )
    # If the current frame count is more than 80% of the stack limit, or comes within 200 frames, increase the
    # stack limit by ``frame_delta``.
    if frame_count > min(0.8 * stack_limit, stack_limit - 200):
        LOGGER.warning(
            "Current stack contains %d frames which is close to the limit of %d. Increasing the limit by %d",
            frame_count,
            stack_limit,
            frame_delta,
        )
        sys.setrecursionlimit(stack_limit + frame_delta)
    manager = get_manager()
    runner = manager.get_runner()
    inputs = process_class.create_inputs(*args, **kwargs)
    # Remove all the known inputs from the kwargs
    for port in process_class.spec().inputs:
        kwargs.pop(port, None)
    # If any kwargs remain, the spec should be dynamic, so we raise if it isn't
    if kwargs and not process_class.spec().inputs.dynamic:
        raise ValueError(
            f"{function.__name__} does not support these kwargs: {kwargs.keys()}"
        )
    if parent_pid:
        process = process_class(inputs=inputs, runner=runner, parent_pid=parent_pid)
    else:
        process = process_class(inputs=inputs, runner=runner)
    # Only add handlers for interrupt signal to kill the process if we are in a local and not a daemon runner.
    # Without this check, running process functions in a daemon worker would be killed if the daemon is shutdown
    current_runner = manager.get_runner()
    original_handler = None
    kill_signal = signal.SIGINT
    if not current_runner.is_daemon_runner:

        def kill_process(_num, _frame):
            """Send the kill signal to the process in the current scope."""
            LOGGER.critical(
                "runner received interrupt, killing process %s", process.pid
            )
            result = process.kill(
                msg="Process was killed because the runner received an interrupt"
            )
            return result

        # Store the current handler on the signal such that it can be restored after process has terminated
        original_handler = signal.getsignal(kill_signal)
        signal.signal(kill_signal, kill_process)
    try:
        result = process.execute()
    finally:
        # If the `original_handler` is set, that means the `kill_process` was bound, which needs to be reset
        if original_handler:
            signal.signal(signal.SIGINT, original_handler)
    store_provenance = inputs.get("metadata", {}).get("store_provenance", True)
    if not store_provenance:
        process.node._storable = False
        process.node._unstorable_message = (
            "cannot store node because it was run with `store_provenance=False`"
        )
    return result, process.node


def instantiate_process(
    runner: "Runner",
    process: Union["Process", Type["Process"], "ProcessBuilder"],
    _parent_pid=None,
    **inputs,
) -> "Process":
    """Return an instance of the process with the given inputs. The function can deal with various types
    of the `process`:

        * Process instance: will simply return the instance
        * ProcessBuilder instance: will instantiate the Process from the class and inputs defined within it
        * Process class: will instantiate with the specified inputs

    If anything else is passed, a ValueError will be raised

    :param process: Process instance or class, CalcJobNode class or ProcessBuilder instance
    :param inputs: the inputs for the process to be instantiated with
    """

    if isinstance(process, Process):
        assert not inputs
        assert runner is process.runner
        return process

    if isinstance(process, ProcessBuilder):
        builder = process
        process_class = builder.process_class
        inputs.update(**builder._inputs(prune=True))
    elif is_process_function(process):
        process_class = process.process_class  # type: ignore[attr-defined]
    elif inspect.isclass(process) and issubclass(process, Process):
        process_class = process
    else:
        raise ValueError(
            f"invalid process {type(process)}, needs to be Process or ProcessBuilder"
        )

    if _parent_pid:
        process = process_class(runner=runner, inputs=inputs, parent_pid=_parent_pid)
    else:
        process = process_class(runner=runner, inputs=inputs)

    return process


def submit(
    process: TYPE_SUBMIT_PROCESS,
    inputs: dict[str, t.Any] | None = None,
    *,
    wait: bool = False,
    wait_interval: int = 5,
    **kwargs: t.Any,
) -> ProcessNode:
    """Submit the process with the supplied inputs to the daemon immediately returning control to the interpreter.

    .. warning: this should not be used within another process. Instead, there one should use the ``submit`` method of
        the wrapping process itself, i.e. use ``self.submit``.

    .. warning: submission of processes requires ``store_provenance=True``.

    :param process: the process class, instance or builder to submit
    :param inputs: the input dictionary to be passed to the process
    :param wait: when set to ``True``, the submission will be blocking and wait for the process to complete at which
        point the function returns the calculation node.
    :param wait_interval: the number of seconds to wait between checking the state of the process when ``wait=True``.
    :param kwargs: inputs to be passed to the process. This is an alternative to the positional ``inputs`` argument.
    :return: the calculation node of the process
    """
    _parent_pid = kwargs.pop("_parent_pid", None)
    runner = kwargs.pop("runner", None)
    inputs = prepare_inputs(inputs, **kwargs)

    # Submitting from within another process requires ``self.submit``` unless it is a work function, in which case the
    # current process in the scope should be an instance of ``FunctionProcess``.
    # if is_process_scoped() and not isinstance(Process.current(), FunctionProcess):
    # raise InvalidOperation('Cannot use top-level `submit` from within another process, use `self.submit` instead')

    if not runner:
        runner = manager.get_manager().get_runner()
    assert runner.persister is not None, "runner does not have a persister"
    assert runner.controller is not None, "runner does not have a controller"

    process_inited = instantiate_process(
        runner, process, _parent_pid=_parent_pid, **inputs
    )

    # If a dry run is requested, simply forward to `run`, because it is not compatible with `submit`. We choose for this
    # instead of raising, because in this way the user does not have to change the launcher when testing. The same goes
    # for if `remote_folder` is present in the inputs, which means we are importing an already completed calculation.
    if process_inited.metadata.get("dry_run", False) or "remote_folder" in inputs:
        _, node = run_get_node(process_inited)
        return node

    if not process_inited.metadata.store_provenance:
        raise InvalidOperation("cannot submit a process with `store_provenance=False`")

    runner.persister.save_checkpoint(process_inited)
    process_inited.close()

    # Do not wait for the future's result, because in the case of a single worker this would cock-block itself
    runner.controller.continue_process(process_inited.pid, nowait=False, no_reply=True)
    node = process_inited.node

    if not wait:
        return node

    while not node.is_terminated:
        LOGGER.report(
            f"Process<{node.pk}> has not yet terminated, current state is `{node.process_state}`. "
            f"Waiting for {wait_interval} seconds."
        )
        time.sleep(wait_interval)

    return node


def start_scheduler_worker(foreground: bool = False) -> None:
    """Start a scheduler worker for the currently configured profile.

    :param foreground: If true, the logging will be configured to write to stdout, otherwise it will be configured to
        write to the scheduler log file.
    """
    import asyncio
    import signal
    import sys

    from aiida.common.log import configure_logging
    from aiida.engine.daemon.client import get_daemon_client
    from aiida.manage import get_config_option, get_manager
    from aiida_workgraph.engine.scheduler import WorkGraphScheduler

    daemon_client = get_daemon_client()
    configure_logging(
        daemon=not foreground, daemon_log_file=daemon_client.daemon_log_file
    )

    LOGGER.debug(f"sys.executable: {sys.executable}")
    LOGGER.debug(f"sys.path: {sys.path}")

    try:
        manager = get_manager()
        # runner = manager.create_daemon_runner()
        runner = manager.create_runner(broker_submit=True)
        manager.set_runner(runner)
    except Exception:
        LOGGER.exception("daemon worker failed to start")
        raise

    if isinstance(rlimit := get_config_option("daemon.recursion_limit"), int):
        LOGGER.info("Setting maximum recursion limit of daemon worker to %s", rlimit)
        sys.setrecursionlimit(rlimit)

    signals = (signal.SIGTERM, signal.SIGINT)
    for s in signals:
        # https://github.com/python/mypy/issues/12557
        runner.loop.add_signal_handler(s, lambda s=s: asyncio.create_task(shutdown_worker(runner)))  # type: ignore[misc]

    print("runner", runner)
    process_inited = instantiate_process(runner, WorkGraphScheduler)
    runner.loop.create_task(process_inited.step_until_terminated())
    print("node", process_inited.node)

    try:
        LOGGER.info("Starting a daemon worker")
        runner.start()
    except SystemError as exception:
        LOGGER.info("Received a SystemError: %s", exception)
        runner.close()

    LOGGER.info("Daemon worker started")
