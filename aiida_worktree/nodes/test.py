from scinode.core.node import Node


class AiiDAInt(Node):
    identifier = "AiiDAInt"
    name = "AiiDAInt"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self):
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self):
        inp = self.inputs.new("General", "value", default=0.0)
        inp.add_property("AiiDAInt", default=1.0)
        self.outputs.new("Int", "result")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "Int",
        }


class AiiDAFloat(Node):
    identifier = "AiiDAFloat"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self):
        self.properties.new("Float", "t", default=1.0)

    def create_sockets(self):
        self.inputs.new("Float", "value", default=0.0)
        self.outputs.new("Float", "result")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "Float",
        }


class AiiDAString(Node):
    identifier = "AiiDAString"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self):
        self.properties.new("Float", "t", default=1.0)

    def create_sockets(self):
        self.inputs.new("String", "value", default="")
        self.outputs.new("String", "result")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "Str",
        }


class AiiDAList(Node):
    identifier = "AiiDAList"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self):
        self.properties.new("BaseList", "value", default=[])

    def create_sockets(self):
        self.outputs.new("General", "Parameters")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "List",
        }


class AiiDADict(Node):
    identifier = "AiiDADict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self):
        self.properties.new("BaseDict", "value", default={})

    def create_sockets(self):
        self.outputs.new("General", "Parameters")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "Dict",
        }


class AiiDANode(Node):
    """AiiDANode"""

    identifier = "AiiDANode"
    name = "AiiDANode"
    node_type = "node"
    catalog = "Test"
    args = ["value"]

    def create_properties(self):
        pass

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "value")
        self.outputs.new("General", "node")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "load_node",
        }


class AiiDACode(Node):
    """AiiDACode"""

    identifier = "AiiDACode"
    name = "AiiDACode"
    node_type = "node"
    catalog = "Test"
    args = ["value"]

    def create_properties(self):
        self.properties.new("General", "value", default=1)

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.outputs.new("General", "Code")

    def get_executor(self):
        return {
            "path": "aiida.orm",
            "name": "load_code",
        }


class AiiDAAdd(Node):

    identifier: str = "AiiDAAdd"
    name = "AiiDAAdd"
    node_type = "calcfunction"
    catalog = "Test"

    args = ["x", "y"]
    kwargs = ["t"]

    def create_properties(self):
        self.properties.new("Float", "t", default=1.0)

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("Float", "x")
        inp.add_property("AiiDAFloat", "x", default=0.0)
        inp = self.inputs.new("Float", "y")
        inp.add_property("AiiDAFloat", "y", default=0.0)
        self.outputs.new("Float", "sum")

    def get_executor(self):
        return {
            "path": "aiida_worktree.executors.test",
            "name": "add",
        }


class AiiDAGreater(Node):

    identifier: str = "AiiDAGreater"
    name = "AiiDAGreater"
    node_type = "calcfunction"
    catalog = "Test"
    kwargs = ["x", "y"]

    def create_properties(self):
        pass

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Float", "x")
        self.inputs.new("Float", "y")
        self.outputs.new("Bool", "result")

    def get_executor(self):
        return {
            "path": "aiida_worktree.executors.test",
            "name": "greater",
        }


class AiiDASumDiff(Node):

    identifier: str = "AiiDASumDiff"
    name = "AiiDASumDiff"
    node_type = "calcfunction"
    catalog = "Test"

    args = ["x", "y"]
    kwargs = ["t"]

    def create_properties(self):
        self.properties.new("Float", "t", default=1.0)

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("Float", "x")
        inp.add_property("AiiDAFloat", "x", default=0.0)
        inp = self.inputs.new("Float", "y")
        inp.add_property("AiiDAFloat", "y", default=0.0)
        self.outputs.new("Float", "sum")
        self.outputs.new("Float", "diff")

    def get_executor(self):
        return {
            "path": "aiida_worktree.executors.test",
            "name": "sum_diff",
        }


class AiiDAArithmeticAdd(Node):

    identifier: str = "AiiDAArithmeticAdd"
    name = "AiiDAArithmeticAdd"
    node_type = "calcjob"
    catalog = "Test"
    kwargs = ["code", "x", "y"]

    def create_properties(self):
        pass

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "code")
        inp = self.inputs.new("Int", "x")
        inp.add_property("AiiDAInt", "x", default=0.0)
        inp = self.inputs.new("Int", "y")
        inp.add_property("AiiDAInt", "y", default=0.0)
        self.outputs.new("Int", "sum")

    def get_executor(self):
        return {
            "path": "aiida.calculations.arithmetic.add",
            "name": "ArithmeticAddCalculation",
        }


class AiiDAArithmeticMultiplyAdd(Node):

    identifier: str = "AiiDAArithmeticMultiplyAdd"
    name = "AiiDAArithmeticMultiplyAdd"
    node_type = "workchain"
    catalog = "Test"
    kwargs = ["code", "x", "y", "z"]

    def create_properties(self):
        pass

    def create_sockets(self):
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "code")
        inp = self.inputs.new("Int", "x")
        inp.add_property("AiiDAInt", "x", default=0.0)
        inp = self.inputs.new("Int", "y")
        inp.add_property("AiiDAInt", "y", default=0.0)
        inp = self.inputs.new("Int", "z")
        inp.add_property("AiiDAInt", "z", default=0.0)
        self.outputs.new("Int", "result")

    def get_executor(self):
        return {
            "name": "core.arithmetic.multiply_add",
            "type": "WorkflowFactory",
        }
