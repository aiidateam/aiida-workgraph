from aiida import orm


type_mapping = {
    "default": "workgraph.any",
    "namespace": "workgraph.namespace",
    int: "workgraph.int",
    float: "workgraph.float",
    str: "workgraph.string",
    bool: "workgraph.bool",
    list: "workgraph.list",
    dict: "workgraph.dict",
    orm.Int: "workgraph.int",
    orm.Float: "workgraph.float",
    orm.Str: "workgraph.string",
    orm.Bool: "workgraph.bool",
    orm.List: "workgraph.list",
    orm.Dict: "workgraph.dict",
    orm.StructureData: "workgraph.aiida_structuredata",
}
