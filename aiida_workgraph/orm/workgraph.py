"""Module with `Node` sub class for work processes."""
from typing import Optional, Tuple

from aiida.common.lang import classproperty

from aiida.orm.nodes.process.workflow.workchain import WorkChainNode

__all__ = ("WorkGraphNode",)


class WorkGraphNode(WorkChainNode):
    """ORM class for all nodes representing the execution of a WorkGraph."""

    WORKTREE_STATE_INFO_KEY = "workgraph_state_info"

    @classproperty
    def _updatable_attributes(cls) -> Tuple[str, ...]:  # type: ignore
        # pylint: disable=no-self-argument
        return super()._updatable_attributes + (cls.WORKTREE_STATE_INFO_KEY,)

    @property
    def workgraph_state_info(self) -> Optional[str]:
        """
        Return the workgraph state info

        :returns: string representation of the workgraph state info
        """
        return self.base.attributes.get(self.WORKTREE_STATE_INFO_KEY, None)

    def set_workgraph_state_info(self, workgraph_state_info: str) -> None:
        """
        Set the workgraph state info

        :param state: string representation of the workgraph state info
        """
        return self.base.attributes.set(
            self.WORKTREE_STATE_INFO_KEY, workgraph_state_info
        )
