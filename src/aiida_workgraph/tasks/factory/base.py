from node_graph.nodes.factory.base import BaseNodeFactory
from aiida_workgraph.task import Task


class BaseTaskFactory(BaseNodeFactory):
    """
    A base factory to create specialized subclasses of Task,
    embedding the 'ndata' (i.e., all relevant data).
    """

    default_base_class = Task
