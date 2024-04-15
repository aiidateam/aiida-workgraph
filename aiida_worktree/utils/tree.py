def node_creation_hook(self, node):
    """Hook for node creation.

    Args:
        node (Node): a node to be created.
    """
    # send message to the widget
    self.parent.send_message("Node created: {}".format(node.name))


def node_deletion_hook(self, node):
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
