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
    AiiDAPW,
    AiiDADos,
    AiiDAProjwfc,
)

node_list = [
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
    AiiDAPW,
    AiiDADos,
    AiiDAProjwfc,
]


# should after node_list, otherwise circular import
node_pool = get_entries(entry_point_name="aiida_worktree.node")
