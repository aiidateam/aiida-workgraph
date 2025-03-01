from .aiida import CalcFunctionTask, WorkChainTask


class TestAdd(CalcFunctionTask):

    identifier: str = "workgraph.test_add"
    name = "TestAAdd"
    node_type = "CALCFUNCTION"
    catalog = "Test"

    def create_properties(self) -> None:
        self.add_property("workgraph.aiida_float", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        inp = self.add_input("workgraph.aiida_float", "x")
        inp.add_property("workgraph.aiida_float", "x", default=0.0)
        inp = self.add_input("workgraph.aiida_float", "y")
        inp.add_property("workgraph.aiida_float", "y", default=0.0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        self.add_output("workgraph.aiida_float", "sum")
        self.add_output("workgraph.any", "_wait")
        self.add_output("workgraph.any", "_outputs")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.test",
            "callable_name": "add",
        }
        return executor


class TestSumDiff(CalcFunctionTask):

    identifier: str = "workgraph.test_sum_diff"
    name = "TestSumDiff"
    catalog = "Test"

    def create_properties(self) -> None:
        self.properties._new("workgraph.aiida_float", "t", default=1.0)

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        inp = self.add_input("workgraph.aiida_float", "x")
        inp.add_property("workgraph.aiida_float", "x", default=0.0)
        inp = self.add_input("workgraph.aiida_float", "y")
        inp.add_property("workgraph.aiida_float", "y", default=0.0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        self.add_output("workgraph.aiida_float", "sum")
        self.add_output("workgraph.aiida_float", "diff")
        self.add_output("workgraph.any", "_wait")
        self.add_output("workgraph.any", "_outputs")

    def get_executor(self):
        executor = {
            "module_path": "aiida_workgraph.executors.test",
            "callable_name": "sum_diff",
        }
        return executor


class TestArithmeticMultiplyAdd(WorkChainTask):

    identifier: str = "workgraph.test_arithmetic_multiply_add"
    name = "TestArithmeticMultiplyAdd"
    catalog = "Test"

    def create_properties(self) -> None:
        pass

    def create_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input("workgraph.any", "code")
        inp = self.add_input("workgraph.aiida_int", "x")
        inp.add_property("workgraph.aiida_int", "x", default=0)
        inp = self.add_input("workgraph.aiida_int", "y")
        inp.add_property("workgraph.aiida_int", "y", default=0)
        inp = self.add_input("workgraph.aiida_int", "z")
        inp.add_property("workgraph.aiida_int", "z", default=0)
        self.add_input(
            "workgraph.any", "_wait", metadata={"arg_type": "none"}, link_limit=100000
        )
        self.add_output("workgraph.aiida_int", "result")
        self.add_output("workgraph.any", "_wait")
        self.add_output("workgraph.any", "_outputs")

    def get_executor(self):
        executor = {
            "module_path": "aiida.workflows.arithmetic.multiply_add",
            "callable_name": "MultiplyAddWorkChain",
        }
        return executor
