from __future__ import annotations
from typing import Callable, List, Optional, Type, Dict, TYPE_CHECKING
from aiida_workgraph.socket_spec import (
    from_aiida_process,
    infer_specs_from_callable,
)
from node_graph.socket_spec import SocketSpec, merge_specs
from node_graph.node_spec import NodeSpec
from node_graph.executor import RuntimeExecutor
from node_graph.error_handler import ErrorHandlerSpec, normalize_error_handlers

if TYPE_CHECKING:
    from node_graph import Node


def build_callable_nodespec(
    *,
    obj: Callable,
    node_type: str,
    base_class: Type['Node'],
    identifier: Optional[str] = None,
    catalog: str = 'Others',
    in_spec: Optional[SocketSpec | List[str]] = None,
    out_spec: Optional[SocketSpec | List[str]] = None,
    process_cls: Optional[type] = None,  # e.g. PythonJob, PyFunction, or aiida.engine.Process
    add_inputs: Optional[SocketSpec | List[str]] = None,
    add_outputs: Optional[SocketSpec | List[str]] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
    metadata: Optional[dict] = None,
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

    # 1) infer from the callable (keep a snapshot before augmentation)
    func_in, func_out = infer_specs_from_callable(obj, in_spec, out_spec)
    # "metadata" is reserved for AiiDA process, so raise error if user tries to use it
    if 'metadata' in func_in.fields:
        fn = getattr(obj, '__name__', 'the task function')
        raise ValueError(
            "Invalid input name: 'metadata'\n"
            "Reason: In AiiDA, 'metadata' is reserved for process-level settings "
            '(e.g., call_link_label, description) and cannot be used as a task input.\n\n'
            f'How to fix: Rename the argument in {fn} to something else, e.g.: task_metadata.\n\n'
            'Example:\n'
            '    # before\n'
            f'    def {fn}(metadata: dict, x: int):\n'
            '        ...\n\n'
            '    # after\n'
            f'    def {fn}(task_metadata: dict, x: int):\n'
            '        ...\n'
        )

    # 2) process-contributed I/O (if any)
    proc_in = proc_out = None
    if process_cls is not None:
        proc_in, proc_out = from_aiida_process(process_cls)
        func_in = merge_specs(func_in, proc_in)
        func_out = merge_specs(func_out, proc_out)

    # 3) additional fields (if any)
    if add_inputs is not None:
        func_in = merge_specs(func_in, add_inputs)
    if add_outputs is not None:
        func_out = merge_specs(func_out, add_outputs)

    # 4) metadata: keep a record of each contribution
    metadata = metadata or {}
    metadata.update(
        {
            'non_function_inputs': list(
                set((proc_in and proc_in.fields.keys()) or []) | set((add_inputs and add_inputs.fields.keys()) or [])
            ),
            'non_function_outputs': list(
                set((proc_out and proc_out.fields.keys()) or [])
                | set((add_outputs and add_outputs.fields.keys()) or [])
            ),
        }
    )

    return NodeSpec(
        identifier=identifier or obj.__name__,
        node_type=node_type,
        catalog=catalog,
        inputs=func_in,
        outputs=func_out,
        executor=RuntimeExecutor.from_callable(obj),
        error_handlers=error_handlers,
        base_class=base_class,
        metadata=metadata,
    )
