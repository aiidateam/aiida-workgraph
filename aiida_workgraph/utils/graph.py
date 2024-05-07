from typing import Any


def node_creation_hook(self, node: Any) -> None:
    """Hook for node creation.

    Args:
        node (Node): a node to be created.
    """
    # send message to the widget
    self.parent._widget.send(
        {"type": "add_node", "data": {"label": node.name, "inputs": [], "outputs": []}}
    )


def node_deletion_hook(self, node: Any) -> None:
    """Hook for node deletion.

    Args:
        node (Node): a node to be deleted.
    """
    # remove all links to the node
    link_index = []
    for index, link in enumerate(self.parent.links):
        if link.from_node.name == node.name or link.to_node.name == node.name:
            link_index.append(index)
    del self.parent.links[link_index]
    self.parent._widget.send({"type": "delete_node", "data": {"name": node.name}})


def link_creation_hook(self, link: Any) -> None:
    """Hook for link creation.

    Args:
        link (Link): a link to be created.
    """
    self.parent._widget.send(
        {
            "type": "add_link",
            "data": {
                "from_node": link.from_node.name,
                "from_socket": link.from_socket.name,
                "to_node": link.to_node.name,
                "to_socket": link.to_socket.name,
            },
        }
    )


def link_deletion_hook(self, link: Any) -> None:
    """Hook for link deletion.

    Args:
        link (Link): a link to be deleted.
    """
    self.parent._widget.send(
        {
            "type": "delete_link",
            "data": {
                "from_node": link.from_node.name,
                "from_socket": link.from_socket.name,
                "to_node": link.to_node.name,
                "to_socket": link.to_socket.name,
            },
        }
    )
