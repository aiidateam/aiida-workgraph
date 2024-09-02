from plumpy.process_comms import RemoteProcessThreadController
from typing import Any, Optional

"""
Note: I modified the the create_daemon_runner function and RemoteProcessThreadController
to include the queue_name argument.

"""


def create_daemon_runner(
    manager, queue_name: str = None, loop: Optional["asyncio.AbstractEventLoop"] = None
) -> "Runner":
    """Create and return a new daemon runner.
    This is used by workers when the daemon is running and in testing.
    :param loop: the (optional) asyncio event loop to use
    :return: a runner configured to work in the daemon configuration
    """
    from plumpy.persistence import LoadSaveContext
    from aiida.engine import persistence
    from aiida.engine.processes.launcher import ProcessLauncher
    from plumpy.communications import convert_to_comm

    runner = manager.create_runner(broker_submit=True, loop=loop)
    runner_loop = runner.loop
    # Listen for incoming launch requests
    task_receiver = ProcessLauncher(
        loop=runner_loop,
        persister=manager.get_persister(),
        load_context=LoadSaveContext(runner=runner),
        loader=persistence.get_object_loader(),
    )

    def callback(_comm, msg):
        print("Received message: {}".format(msg))
        import asyncio

        asyncio.run(task_receiver(_comm, msg))
        print("task_receiver._continue done")
        return True

    assert runner.communicator is not None, "communicator not set for runner"
    if queue_name is not None:
        print("queue_name: {}".format(queue_name))
        queue = runner.communicator._communicator.task_queue(
            queue_name, prefetch_count=1
        )
        # queue.add_task_subscriber(callback)
        # important to convert the callback
        converted = convert_to_comm(task_receiver, runner.communicator._loop)
        queue.add_task_subscriber(converted)
    else:
        runner.communicator.add_task_subscriber(task_receiver)
    return runner


class ControllerWithQueueName(RemoteProcessThreadController):
    def __init__(self, queue_name: str, **kwargs):
        super().__init__(**kwargs)
        self.queue_name = queue_name

    def task_send(self, message: Any, no_reply: bool = False) -> Optional[Any]:
        """
        Send a task to be performed using the communicator

        :param message: the task message
        :param no_reply: if True, this call will be fire-and-forget, i.e. no return value
        :return: the response from the remote side (if no_reply=False)
        """
        queue = self._communicator.task_queue(self.queue_name)
        return queue.task_send(message, no_reply=no_reply)
