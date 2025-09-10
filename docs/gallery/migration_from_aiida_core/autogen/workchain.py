"""
=================================
Convert WorkChain to WorkGraph
=================================

This tutorial demonstrates how to convert an AiiDA WorkChain into a WorkGraph.
We start by defining a simple WorkChain that sums even numbers and then
re-implement its logic as a WorkGraph.

"""
# %%
# WorkChain
# ===========================
# First, we define a standard AiiDA WorkChain. This example, `SumEvenWorkChain`,
# calculates the sum of all even integers from 1 up to a given number `N`.
# It uses a `while` loop to iterate and an `if` condition to check for even numbers.

from aiida.engine import WorkChain, calcfunction, if_, while_
from aiida.orm import Int


@calcfunction
def add(x: Int, y: Int) -> Int:
    """A simple calcfunction to add two integers."""
    return Int(x.value + y.value)


class SumEvenWorkChain(WorkChain):
    """WorkChain to sum all even numbers from 1 up to N."""

    @classmethod
    def define(cls, spec):
        """Specify inputs, outputs, and the workchain logic."""
        super().define(spec)
        spec.input("N", valid_type=Int, help="The integer to sum up to.")
        spec.outline(
            cls.setup,
            while_(cls.smaller_than)(
                if_(cls.is_even)(
                    cls.add_total,
                ),
                cls.update_n,
            ),
            cls.result,
        )
        spec.output("total", valid_type=Int, help="The final sum.")

    def setup(self):
        """Initialize context variables."""
        self.ctx.n = Int(1)
        self.ctx.total = Int(0)

    def smaller_than(self):
        """Condition for the while loop: check if n < N."""
        return self.ctx.n.value < self.inputs.N.value

    def is_even(self):
        """Condition for the if statement: check if n is even."""
        return self.ctx.n.value % 2 == 0

    def add_total(self):
        """Add the current number to the total."""
        self.ctx.total = add(self.ctx.total, self.ctx.n)

    def update_n(self):
        """Increment the current number."""
        self.ctx.n = add(self.ctx.n, Int(1))

    def result(self):
        """Attach the final sum to the outputs."""
        self.out("total", self.ctx.total)


# %%
# WorkGraph Equivalent
# ========================
# Now, we convert the `SumEvenWorkChain` into a `WorkGraph`. The core logic
# remains the same, but the implementation differs. The outline of the
# `SumEvenWorkChain` provides a clear blueprint for our graph:
#
# .. code-block:: python
#
#    spec.outline(
#        cls.setup,
#        while_(cls.smaller_than)(
#            if_(cls.is_even)(
#                cls.add_total,
#            ),
#            cls.update_n,
#        ),
#        cls.result,
#        )
#
# Each method in the outline is transformed into a `task`.

from aiida_workgraph import task, While, If

# First, we convert the existing `add` calcfunction into a reusable task.
add_task = task(add)

# Next, we create tasks for the conditions in our loops and conditionals.
@task()
def smaller_than(n: int, N: int) -> bool:
    """Task to check if n < N."""
    return n < N


@task()
def is_even(n: int) -> bool:
    """Task to check if n is even."""
    return n % 2 == 0


# Finally, we define the WorkGraph itself.
@task.graph()
def sum_even_workgraph(N: int):
    """WorkGraph to sum all even numbers from 1 up to N."""
    from aiida_workgraph.manager import get_current_graph

    wg = get_current_graph()

    # The 'setup' step: initialize context variables.
    wg.ctx = {"n": 1, "total": 0}

    # The 'while' loop. The condition is now a task.
    with While(smaller_than(wg.ctx.n, N).result):
        # The 'if' condition.
        with If(is_even(wg.ctx.n).result) as if_zone:
            # The 'add_total' step.
            wg.ctx.total = add_task(wg.ctx.total, wg.ctx.n).result

        # The 'update_n' step.
        n_new = add_task(wg.ctx.n, 1)

        # Manually set a dependency to ensure the 'if' block completes
        # before 'n' is updated for the next iteration.
        if_zone >> n_new
        wg.ctx.n = n_new.result

    # The 'result' step: define the final output of the graph.
    return wg.ctx.total


# %%
# Running the WorkGraph
# =====================
# With the WorkGraph defined, we can now generate it, inspect its structure,
# and execute it to get the result.

from aiida import load_profile

# Load your AiiDA profile.
load_profile()

# Generate the WorkGraph instance with a specific input.
N = 5
wg = sum_even_workgraph.build(N=N)

# The `to_html()` method generates an interactive visualization of the graph.
# In a Sphinx-Gallery build, this will be embedded directly in the output.
wg.to_html()

# %%
# Execute the WorkGraph and print the result.
wg.run()
print(f"The sum of even numbers up to {N} is: {wg.outputs['result'].value}")


# %%
# Conclusion
# ==========
# This tutorial has shown the process of converting an AiiDA
# WorkChain to a WorkGraph, by mapping the procedural steps of a WorkChain
# (setup, loop, condition, action, result) to a graph of interconnected tasks.
