from node_graph.collection import EntryPointPool

# global instance
SocketPool = EntryPointPool(entry_point_group="aiida_workgraph.socket")
SocketPool["any"] = SocketPool.workgraph.any
SocketPool["namespace"] = SocketPool.workgraph.namespace
