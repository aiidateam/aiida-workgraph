from typing import Dict
from aiida_workgraph.task import Task


class AiiDAInt(Task):
    identifier = "AiiDAInt"
    name = "AiiDAInt"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self) -> None:
        inp = self.inputs.new("Any", "value", default=0.0)
        inp.add_property("AiiDAInt", default=1.0)
        self.outputs.new("AiiDAInt", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Int",
        }


class AiiDAFloat(Task):
    identifier = "AiiDAFloat"
    name = "AiiDAFloat"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs.new("AiiDAFloat", "value", default=0.0)
        self.outputs.new("AiiDAFloat", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Float",
        }


class AiiDAString(Task):
    identifier = "AiiDAString"
    name = "AiiDAString"
    node_type = "data"
    catalog = "Test"

    args = ["value"]
    kwargs = ["t"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs.new("AiiDAString", "value", default="")
        self.outputs.new("AiiDAString", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Str",
        }


class AiiDAList(Task):
    identifier = "AiiDAList"
    name = "AiiDAList"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseList", "value", default=[])

    def create_sockets(self) -> None:
        self.outputs.new("Any", "Parameters")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "List",
        }


class AiiDADict(Task):
    identifier = "AiiDADict"
    name = "AiiDADict"
    node_type = "data"
    catalog = "Test"

    args = ["value"]

    def create_properties(self) -> None:
        self.properties.new("BaseDict", "value", default={})

    def create_sockets(self) -> None:
        self.outputs.new("Any", "Parameters")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "Dict",
        }


class AiiDANode(Task):
    """AiiDANode"""

    identifier = "AiiDANode"
    name = "AiiDANode"
    node_type = "node"
    catalog = "Test"
    kwargs = ["identifier", "pk", "uuid", "label"]

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Any", "identifier")
        self.inputs.new("Any", "pk")
        self.inputs.new("Any", "uuid")
        self.inputs.new("Any", "label")
        self.outputs.new("Any", "node")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "load_node",
        }


class AiiDACode(Task):
    """AiiDACode"""

    identifier = "AiiDACode"
    name = "AiiDACode"
    node_type = "node"
    catalog = "Test"
    kwargs = ["identifier", "pk", "uuid", "label"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Any", "identifier")
        self.inputs.new("Any", "pk")
        self.inputs.new("Any", "uuid")
        self.inputs.new("Any", "label")
        self.outputs.new("Any", "Code")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "load_code",
        }


class AiiDAAdd(Task):

    identifier: str = "AiiDAAdd"
    name = "AiiDAAdd"
    node_type = "CALCFUNCTION"
    catalog = "Test"

    args = ["x", "y"]
    kwargs = ["t"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("AiiDAFloat", "x")
        inp.add_property("AiiDAFloat", "x", default=0.0)
        inp = self.inputs.new("AiiDAFloat", "y")
        inp.add_property("AiiDAFloat", "y", default=0.0)
        self.outputs.new("AiiDAFloat", "sum")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.test",
            "name": "add",
        }


class AiiDAGreater(Task):

    identifier: str = "AiiDAGreater"
    name = "AiiDAGreater"
    node_type = "CALCFUNCTION"
    catalog = "Test"
    kwargs = ["x", "y"]

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("AiiDAFloat", "x")
        self.inputs.new("AiiDAFloat", "y")
        self.outputs.new("AiiDABool", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.test",
            "name": "greater",
        }


class AiiDASumDiff(Task):

    identifier: str = "AiiDASumDiff"
    name = "AiiDASumDiff"
    node_type = "CALCFUNCTION"
    catalog = "Test"

    args = ["x", "y"]
    kwargs = ["t"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAFloat", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("AiiDAFloat", "x")
        inp.add_property("AiiDAFloat", "x", default=0.0)
        inp = self.inputs.new("AiiDAFloat", "y")
        inp.add_property("AiiDAFloat", "y", default=0.0)
        self.outputs.new("AiiDAFloat", "sum")
        self.outputs.new("AiiDAFloat", "diff")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.test",
            "name": "sum_diff",
        }


class AiiDAArithmeticMultiplyAdd(Task):

    identifier: str = "AiiDAArithmeticMultiplyAdd"
    name = "AiiDAArithmeticMultiplyAdd"
    node_type = "WORKCHAIN"
    catalog = "Test"
    kwargs = ["code", "x", "y", "z"]

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("Any", "code")
        inp = self.inputs.new("AiiDAInt", "x")
        inp.add_property("AiiDAInt", "x", default=0.0)
        inp = self.inputs.new("AiiDAInt", "y")
        inp.add_property("AiiDAInt", "y", default=0.0)
        inp = self.inputs.new("AiiDAInt", "z")
        inp.add_property("AiiDAInt", "z", default=0.0)
        self.outputs.new("AiiDAInt", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "name": "core.arithmetic.multiply_add",
            "type": "WorkflowFactory",
        }
