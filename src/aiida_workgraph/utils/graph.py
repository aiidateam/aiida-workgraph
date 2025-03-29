from typing import Any
from node_graph.collection import NodeCollection


def task_creation_hook(self: NodeCollection, task: Any) -> None:
    """Hook for task creation.

    Args:
        task (Task): a task to be created.
    """
    if not self.parent.interactive_widget:
        return
    # send message to the widget
    if self.parent.widget is not None:
        self.parent.widget.send(
            {
                "type": "add_task",
                "data": {"label": task.name, "inputs": [], "outputs": []},
            }
        )


def task_deletion_hook(self: NodeCollection, task: Any) -> None:
    """Hook for task deletion.

    Args:
        task (Task): a task to be deleted.
    """
    if not self.parent.interactive_widget:
        return
    # remove all links to the task
    link_index = []
    for index, link in enumerate(self.parent.links):
        if link.from_node.name == task.name or link.to_node.name == task.name:
            link_index.append(index)
    del self.parent.links[link_index]
    if self.widget is not None:
        self.parent.widget.send({"type": "delete_node", "data": {"name": task.name}})


def link_creation_hook(self: NodeCollection, link: Any) -> None:
    """Hook for link creation.

    Args:
        link (Link): a link to be created.
    """
    if not self.parent.interactive_widget:
        return
    if self.parent.widget is not None:
        self.parent.widget.send(
            {
                "type": "add_link",
                "data": {
                    "from_node": link.from_node.name,
                    "from_socket": link.from_socket._name,
                    "to_node": link.to_node.name,
                    "to_socket": link.to_socket._name,
                },
            }
        )


def link_deletion_hook(self: NodeCollection, link: Any) -> None:
    """Hook for link deletion.

    Args:
        link (Link): a link to be deleted.
    """
    if not self.parent.interactive_widget:
        return
    if self.parent.widget is not None:
        self.parent.widget.send(
            {
                "type": "delete_link",
                "data": {
                    "from_node": link.from_node.name,
                    "from_socket": link.from_socket._name,
                    "to_node": link.to_node.name,
                    "to_socket": link.to_socket._name,
                },
            }
        )
