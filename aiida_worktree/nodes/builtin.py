from aiida_worktree.node import Node


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


class AiiDAToCtx(Node):
    """AiiDAToCtx"""

    identifier = "ToCtx"
    name = "ToCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key", "value"]

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "key")
        self.inputs.new("General", "value")
        self.outputs.new("General", "result")

    def get_executor(self):
        return {
            "path": "builtins",
            "name": "setattr",
        }


class AiiDAFromCtx(Node):
    """AiiDAFromCtx"""

    identifier = "FromCtx"
    name = "FromCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key"]

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "key")
        self.outputs.new("General", "result")

    def get_executor(self):
        return {
            "path": "builtins",
            "name": "getattr",
        }


if __name__ == "__main__":
    print()
