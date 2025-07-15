from node_graph.collection import EntryPointPool

# global instance
TaskPool = EntryPointPool(entry_point_group="aiida_workgraph.task")
TaskPool["any"] = TaskPool.workgraph.task
TaskPool["graph_inputs"] = TaskPool.workgraph.graph_inputs
TaskPool["graph_outputs"] = TaskPool.workgraph.graph_outputs
TaskPool["graph_ctx"] = TaskPool.workgraph.graph_ctx
