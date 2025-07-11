"""
============================
Migrating from AiiDA Core
============================

"""

# %%
# Introduction
# ============
# In this tutorial, we detail how to convert a WorkChain to a WorkGraph.
# First, we create a example WorkChain in AiiDA Core, then we convert it to a WorkGraph.

from aiida.engine import WorkChain, calcfunction, if_, while_
from aiida.orm import Int
from aiida_workgraph.utils import generate_node_graph


@calcfunction
def add(x, y):
    return Int(x + y)


@calcfunction
def multiply(x, y):
    return Int(x * y)


class AddMultiplyWorkChain(WorkChain):
    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input("x")
        spec.input("y")
        spec.input("z")
        spec.outline(
            cls.setup,
            if_(should_run_add)(
                cls.add,
                cls.inspect_add,
            ),
            while_(smaller_than_n)(
                cls.multiply,
            ),
            cls.results,
        )
        spec.output("result")

    def add(self):
        self.ctx.sum = add(self.inputs.x, self.inputs.y)

    def multiply(self):
        self.ctx.product = multiply(self.ctx.sum, self.inputs.z)

    def results(self):
        self.out("result", self.ctx.product)
