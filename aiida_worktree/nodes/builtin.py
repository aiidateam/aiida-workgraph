from scinode.core.node import Node
from aiida_worktree.executors.builtin import GatherWorkChain


class AiiDAGather(Node):
    """AiiDAGather"""

    identifier = "AiiDAGather"
    name = "AiiDAGather"
    node_type = "workchain"
    catalog = "AiiDA"
    kwargs = ["datas"]

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("General", "datas")
        inp.link_limit = 100000
        self.outputs.new("General", "result")

    def get_executor(self):
        return {
            "path": "aiida_worktree.executors.builtin",
            "name": "GatherWorkChain",
        }


if __name__ == "__main__":
    print(gather_node)
