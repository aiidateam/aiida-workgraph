from node_graph.collection import EntryPointPool

# global instance
PropertyPool = EntryPointPool(entry_point_group="aiida_workgraph.property")
PropertyPool["any"] = PropertyPool.workgraph.any
