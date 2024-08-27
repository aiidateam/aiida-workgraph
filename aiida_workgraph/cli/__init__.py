"""Sub commands of the ``verdi`` command line interface.

The commands need to be imported here for them to be registered with the top-level command group.
"""
from aiida_workgraph.cli import cmd_graph
from aiida_workgraph.cli import cmd_web
from aiida_workgraph.cli import cmd_task
from aiida_workgraph.cli import cmd_scheduler


__all__ = ["cmd_graph", "cmd_web", "cmd_task", "cmd_scheduler"]
