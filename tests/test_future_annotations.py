"""Regression test for #783: dynamic-namespace tasks under PEP 563.

The module-level ``from __future__ import annotations`` is what makes annotations
strings, the condition that breaks run-time signature re-inference of pickled task
callables.
"""

from __future__ import annotations

from typing import Annotated

from aiida_workgraph import dynamic, namespace, task


def test_dynamic_namespace_graph_task_runs_under_future_annotations():
    """A dynamic-namespace ``@task.graph`` must run under PEP 563.

    The engine re-infers the graph task's signature from the unpickled callable;
    under PEP 563 cloudpickle drops the annotation's names, which used to raise
    ``NameError: name 'Annotated' is not defined`` mid-run (#783).
    """

    @task
    def source() -> Annotated[dict, namespace(data=dynamic(int))]:
        return {'data': {'k1': 1, 'k2': 2, 'k3': 3}}

    @task
    def total(data: Annotated[dict, dynamic(int)]) -> int:
        return sum(data.values())

    @task.graph
    def fan(data: Annotated[dict, dynamic(int)]) -> int:
        return total(data=data).result

    @task.graph
    def top() -> int:
        return fan(data=source().data).result

    wg = top.build()
    wg.run()

    assert wg.state == 'FINISHED'
    assert wg.outputs.result.value == 6
