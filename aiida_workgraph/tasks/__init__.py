from node_graph.utils import get_entries

# should after task_list, otherwise circular import
task_pool = {
    **get_entries(entry_point_name="aiida_workgraph.task"),
    **get_entries(entry_point_name="node_graph.task"),
}
