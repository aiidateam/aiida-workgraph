from node_graph.collection import EntryPointPool

# global instance
TaskPool = EntryPointPool(entry_point_group="aiida_workgraph.task")
