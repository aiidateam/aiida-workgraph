from __future__ import annotations
from dataclasses import replace
from copy import deepcopy
from typing import Callable, List, Optional, Any, Type, Dict
from aiida_workgraph.socket_spec import (
    from_aiida_process,
    infer_specs_from_callable,
    namespace,
)
from node_graph.socket_spec import SocketSpec, merge_spec
from node_graph.node_spec import NodeSpec
from node_graph.executor import NodeExecutor
from node_graph.error_handler import ErrorHandlerSpec, normalize_error_handlers

# --- tiny utilities ----------------------------------------------------------


def namespace_with_defaults(defaults: dict[str, Any], **fields: Any) -> SocketSpec:
    """Build a namespace and attach per-field default values."""
    ns = namespace(**fields)
    return replace(ns, defaults=deepcopy(defaults))


def _record_specs_block(
    title: str, in_spec: SocketSpec | None, out_spec: SocketSpec | None
) -> dict:
    """Serialize specs into a metadata block (omit when None)."""
    if in_spec is None and out_spec is None:
        return {}
    block = {}
    if in_spec is not None:
        block["inputs"] = in_spec.to_dict()
    if out_spec is not None:
        block["outputs"] = out_spec.to_dict()
    return {title: block}


# --- one generic builder for callable-based nodes ----------------------------


def build_callable_nodespec(
    *,
    obj: Callable,
    node_type: str,
    base_class: Type["Node"],
    identifier: Optional[str] = None,
    catalog: str = "AIIDA",
    in_spec: Optional[SocketSpec | List[str]] = None,
    out_spec: Optional[SocketSpec | List[str]] = None,
    process_cls: Optional[
        type
    ] = None,  # e.g. PythonJob, PyFunction, or aiida.engine.Process
    add_inputs: Optional[SocketSpec | List[str]] = None,
    add_outputs: Optional[SocketSpec | List[str]] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    executor_factory=None,  # defaults to NodeExecutor.from_callable
) -> NodeSpec:
    """
    - infers function I/O
    - optionally merges process-contributed I/O
    - optionally merges additional I/O
    - records *each* contribution in metadata
    """
    from aiida_workgraph.socket_spec import validate_socket_data

    error_handlers = normalize_error_handlers(error_handlers)

    in_spec = validate_socket_data(in_spec)
    out_spec = validate_socket_data(out_spec)
    executor_factory = executor_factory or NodeExecutor.from_callable

    # 1) infer from the callable (keep a snapshot before augmentation)
    func_in, func_out = infer_specs_from_callable(obj, in_spec, out_spec)
    origin_in, origin_out = deepcopy(func_in), deepcopy(func_out)

    # 2) process-contributed I/O (if any)
    proc_in = proc_out = None
    if process_cls is not None:
        proc_in, proc_out = from_aiida_process(process_cls)
        func_in = merge_spec(func_in, proc_in)
        func_out = merge_spec(func_out, proc_out)

    # 3) additional fields (if any)
    if add_inputs is not None:
        func_in = merge_spec(func_in, add_inputs)
    if add_outputs is not None:
        func_out = merge_spec(func_out, add_outputs)

    # 4) metadata: keep a full audit trail of contributions
    meta_specs = {}
    meta_specs.update(_record_specs_block("function", origin_in, origin_out))
    meta_specs.update(_record_specs_block("process", proc_in, proc_out))
    meta_specs.update(_record_specs_block("additional", add_inputs, add_outputs))

    meta = {
        "node_type": node_type,
        "specs": meta_specs,
    }

    return NodeSpec(
        identifier=identifier or obj.__name__,
        catalog=catalog,
        inputs=func_in,
        outputs=func_out,
        executor=executor_factory(obj),
        error_handlers=error_handlers,
        base_class=base_class,
        metadata=meta,
    )
