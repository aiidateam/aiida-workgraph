from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional


@dataclass(frozen=True)
class MissingInput:
    """One required input socket that has neither a value nor an incoming link.

    Attributes:
        socket_path (str): Dotted path of the socket, e.g. ``'add1.x'`` or ``'graph_inputs.structure'``.
        identifier (str): Socket type identifier, e.g. ``'workgraph.int'`` or ``'workgraph.any'``.
        help (Optional[str]): Help text declared for the socket, if any.
    """

    socket_path: str
    identifier: str
    help: Optional[str] = None

    def __str__(self) -> str:
        return self.socket_path


class MissingRequiredInputsError(ValueError):
    """Raised when a WorkGraph is run with required inputs unset.

    The human-readable message lists the missing sockets; the ``missing`` attribute
    carries them as structured :class:`MissingInput` entries so callers can act on
    them programmatically (e.g. point users to the right input by type).
    """

    def __reduce__(self) -> tuple[type[MissingRequiredInputsError], tuple[List[MissingInput]]]:
        """Reduce to the pattern that ``__init__`` can digest."""
        return (self.__class__, (self.missing,))

    def __init__(self, missing: List[MissingInput]) -> None:
        self.missing = sorted(missing, key=lambda entry: entry.socket_path)
        bullets = '\n'.join(f'  • {entry.socket_path}' for entry in self.missing)
        super().__init__(
            'Missing required inputs:\n'
            f'{bullets}\n\n'
            'How to fix:\n'
            '  1) Provide these values (at build time or by linking from upstream task outputs).\n'
            '  2) If some are intentionally unused, exclude them from the namespace at the call site, e.g.:\n'
            '     Annotated[dict, some_task.inputs, SocketSpecSelect(exclude=["pw.structure", ...])]\n\n'
            "Note: exclude paths are relative to the task's input namespace (e.g. 'pw.structure')."
        )
