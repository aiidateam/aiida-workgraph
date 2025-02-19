from aiida_workgraph.socket import TaskSocket


def if_(condition):
    """Helper function to create an If object."""
    return If_(condition)


class If_:
    def __init__(self, condition: TaskSocket):
        self.condition = condition
        self.wg = condition._parent._node.parent
        self.true_zone = self.wg.add_task("If", name="if_true", conditions=condition)
        self.false_zone = None

    def __call__(self, *tasks):
        """
        Called when the user does:
            if_(condition)(<tasks-for-true>)
        """
        self.true_zone.children.add([*tasks])
        return self

    def else_(self, *tasks):
        """
        Called when the user does:
            .else_(<tasks-for-false>)
        """
        self.false_zone = self.wg.add_task(
            "If", name="if_false", conditions=self.condition, invert_condition=True
        )
        self.false_zone.children.add([*tasks])

        return self
