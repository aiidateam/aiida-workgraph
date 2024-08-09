from node_graph.utils import get_entries

socket_pool = {
    **get_entries(entry_point_name="node_graph.socket"),
    **get_entries(entry_point_name="aiida_workgraph.socket"),
}
