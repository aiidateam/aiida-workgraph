from aiida import orm


type_mapping = {
    "default": "workgraph.any",
    "namespace": "workgraph.namespace",
    int: "workgraph.int",
    float: "workgraph.float",
    str: "workgraph.string",
    bool: "workgraph.bool",
    orm.Int: "workgraph.aiida_int",
    orm.Float: "workgraph.aiida_float",
    orm.Str: "workgraph.aiida_string",
    orm.Bool: "workgraph.aiida_bool",
    orm.List: "workgraph.aiida_list",
    orm.Dict: "workgraph.aiida_dict",
    orm.StructureData: "workgraph.aiida_structuredata",
}
