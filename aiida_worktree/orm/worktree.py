"""Module with `Node` sub class for work processes."""
from typing import Optional, Tuple

from aiida.common.lang import classproperty

from aiida.orm.nodes.process.workflow.workflow import WorkflowNode

__all__ = ("WorkTreeNode",)


class WorkTreeNode(WorkflowNode):
    """ORM class for all nodes representing the execution of a WorkTree."""

    WORKTREE_STATE_INFO_KEY = "worktree_state_info"

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # type: ignore
        # pylint: disable=no-self-argument
        return super()._updatable_attributes + (cls.WORKTREE_STATE_INFO_KEY,)

    @property
    def worktree_state_info(self) -> Optional[str]:
        """
        Return the worktree state info

        :returns: string representation of the worktree state info
        """
        return self.base.attributes.get(self.WORKTREE_STATE_INFO_KEY, None)

    def set_worktree_state_info(self, worktree_state_info: str) -> None:
        """
        Set the worktree state info

        :param state: string representation of the worktree state info
        """
        return self.base.attributes.set(
            self.WORKTREE_STATE_INFO_KEY, worktree_state_info
        )
