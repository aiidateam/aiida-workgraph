from node_graph.utils import get_entries
from .builtin import AiiDAGather, AiiDAToCtx, AiiDAFromCtx
from .test import (
    AiiDAInt,
    AiiDAFloat,
    AiiDAString,
    AiiDAList,
    AiiDADict,
    AiiDANode,
    AiiDACode,
    AiiDAAdd,
    AiiDAGreater,
    AiiDASumDiff,
    AiiDAArithmeticMultiplyAdd,
)
from .qe import (
    AiiDAKpoint,
    AiiDAPWPseudo,
    AiiDAStructure,
)

task_list = [
    AiiDAGather,
    AiiDAToCtx,
    AiiDAFromCtx,
    AiiDAInt,
    AiiDAFloat,
    AiiDAString,
    AiiDAList,
    AiiDADict,
    AiiDANode,
    AiiDACode,
    AiiDAAdd,
    AiiDAGreater,
    AiiDASumDiff,
    AiiDAArithmeticMultiplyAdd,
    AiiDAKpoint,
    AiiDAPWPseudo,
    AiiDAStructure,
]


# should after task_list, otherwise circular import
task_pool = get_entries(entry_point_name="aiida_workgraph.task")
