from aiida.engine import ToContext, WorkChain
from aiida.workflows.arithmetic.multiply_add import MultiplyAddWorkChain
from aiida.calculations.arithmetic.add import ArithmeticAddCalculation
from aiida.common import AttributeDict
from aiida.orm import Int


class WorkChainWithNestNamespace(WorkChain):
    """WorkChain to add two numbers."""

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input_namespace("non_dynamic_port")
        spec.input("non_dynamic_port.a", valid_type=Int)
        spec.input("non_dynamic_port.b", valid_type=Int)
        spec.input("non_dynamic_port.a")
        spec.expose_inputs(
            ArithmeticAddCalculation,
            namespace="add",
        )
        spec.expose_inputs(MultiplyAddWorkChain, namespace="multiply_add")
        spec.outline(
            cls.add,
            cls.multiply_add,
            cls.validate_result,
            cls.result,
        )
        spec.output("result", valid_type=Int)
        spec.expose_outputs(MultiplyAddWorkChain, namespace="multiply_add")
        spec.exit_code(
            400, "ERROR_NEGATIVE_NUMBER", message="The result is a negative number."
        )

    def add(self):
        """Add two numbers using the `ArithmeticAddCalculation` calculation job plugin."""
        inputs = AttributeDict(self.exposed_inputs(ArithmeticAddCalculation, "add"))
        future = self.submit(ArithmeticAddCalculation, **inputs)
        self.report(f"Submitted the `ArithmeticAddCalculation`: {future}")
        return ToContext(addition=future)

    def multiply_add(self):
        """Multiply and add two numbers using the `MultiplyAddWorkChain` workchain."""
        inputs = self.exposed_inputs(MultiplyAddWorkChain, "multiply_add")
        inputs["z"] = self.ctx.addition.outputs.sum
        future = self.submit(MultiplyAddWorkChain, **inputs)
        self.report(f"Submitted the `MultiplyAddWorkChain`: {future}")
        return ToContext(multiply_add=future)

    def validate_result(self):
        """Make sure the result is not negative."""
        result = self.ctx.addition.outputs.sum
        if result.value < 0:
            return self.exit_codes.ERROR_NEGATIVE_NUMBER

    def result(self):
        """Add the result to the outputs."""
        self.out("result", self.ctx.addition.outputs.sum)


class WorkChainWithDynamicNamespace(WorkChain):
    """WorkChain with dynamic namespace."""

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input_namespace("dynamic_port", dynamic=True)
