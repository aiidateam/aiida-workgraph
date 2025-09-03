from node_graph.registry import RegistryHub


registry_hub = RegistryHub.from_prefix(
    node_group="aiida_workgraph.task",
    socket_group="aiida_workgraph.socket",
    property_group="aiida_workgraph.property",
    type_mapping_group="aiida_workgraph.type_mapping",
    identifier_prefix="workgraph",
)
